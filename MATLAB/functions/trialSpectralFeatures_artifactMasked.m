function [trialSpectFeatures, STFTs] = trialSpectralFeatures(data, stim, srate, window, noverlap, nfft, bands, baselinePSD)
% trialSpectralFeatures Extract trial-wise spectral features from multichannel data.
% Syntax
%   trialFeatures = trialSpectralFeatures(data, stim, srate)
%   trialFeatures = trialSpectralFeatures(data, stim, srate, window, noverlap, freqVec, bands, baselinePSD)
%
% Inputs
%   data        - [samples x channels] signal matrix.
%   stim        - [samples x 1] stimulus/state vector.
%                 Assumption: stim == 0 means no active trial and
%                 stim == -1 marks artifact samples that are excluded.
%                 Changes between non-artifact stimulus codes define trials.
%   srate       - Sampling rate in Hz.
% Optional inputs
%   window      - Welch/spectral window. Default: Hann window of ~1 s.
%   noverlap    - Overlap in samples. Default: 50% window overlap.
%   nfft        - number of FFT points. Default = window length 
%                 Default: 1:srate/2.
%   bands       - Struct defining frequency bands.
%                 Example: bands.beta = [15 30].
%   baselinePSD - Optional baseline PSD for normalization. If empty, normalization uses whole-experiment mean PSD.
%                 Accepted dimensions:
%                   [freq x channels]               
% Output
%   trialFeatures - struct with fields:
%       .trialPSD             [channels x frequencies x trials]
%       .meanPSD              [channels x frequencies]
%       .trialNormPSD         [channels x frequencies x trials]
%       .TFmap                [frequencies x time x channels]
%       .TFtime               spectrogram time axis
%       .stimCode             [trials x 1]
%       .validTrial           [trials x 1] logical
%       .trialStartIdx        [trials x 1]
%       .trialStopIdx         [trials x 1]
%       .trialStartTime       [trials x 1] seconds
%       .trialStopTime        [trials x 1] seconds
%       .freqBins             frequency vector
%       .bandPower            struct with [channels x trials] per band
%       .trialCoherence       cell array, one [ch x ch x freq] matrix per trial
%       .info                 parameters used for computation
%
% Notes
%   - PSDs are returned in PSD units, e.g. µV^2/Hz if input is µV.
%   - Band power is integrated over frequency and therefore has units µV^2.
%   - Artifact samples (stim == -1) are excluded before trial-wise analysis.
%   - Trials shorter than the window length after artifact removal are marked invalid.

fprintf('Calculating trial spectral features...\n');

%% Defaults
if nargin < 4 || isempty(window)
    window = hann(2^nextpow2(srate), 'periodic');
end
if nargin < 5 || isempty(noverlap)
    noverlap = floor(length(window)/2);
end
if nargin < 6 || isempty(nfft)
    nfft = length(window);
end
if nargin < 7 || isempty(bands)
    bands = struct( ...
        'delta', [1 4], ...
        'theta', [4 8], ...
        'alpha', [8 15], ...
        'beta', [15 30], ...
        'broadband', [70 115]);
end
if nargin < 8
    baselinePSD = [];
end

%% Validation
data = double(data);
stim = double(stim(:));

if size(data,1) ~= numel(stim)
    error('data and stim must have the same number of samples.');
end

if noverlap >= length(window)
    error('noverlap must be smaller than length(window).');
end

validateBands(bands);
nCh = size(data,2);

%% Catalogue trials while excluding artifact periods (stim == -1)
artifactMask = stim == -1;
validStimIdx = find(~artifactMask);

if isempty(validStimIdx)
    error('No non-artifact samples remain after excluding stim == -1.');
end

trial_sc_vec = zeros(size(stim));      % sample-wise trial number; artifacts remain 0
tr_sc = [];                            % one stim code per trial
trialStartIdx = [];
trialStopIdx = [];

trialCounter = 1;
firstIdx = validStimIdx(1);
trial_sc_vec(firstIdx) = trialCounter;
tr_sc(trialCounter,1) = stim(firstIdx);
trialStartIdx(trialCounter,1) = firstIdx;
prevIdx = firstIdx;
prevStim = stim(firstIdx);

for k = 2:numel(validStimIdx)
    n = validStimIdx(k);

    % Start a new trial when the non-artifact stimulus code changes.
    % Artifact samples between equal stimulus codes are simply clipped out.
    if stim(n) ~= prevStim
        trialStopIdx(trialCounter,1) = prevIdx;
        trialCounter = trialCounter + 1;
        trialStartIdx(trialCounter,1) = n;
        tr_sc(trialCounter,1) = stim(n);
    end

    trial_sc_vec(n) = trialCounter;
    prevIdx = n;
    prevStim = stim(n);
end

% Close final trial
trialStopIdx(trialCounter,1) = prevIdx;
nTrials = trialCounter;
trialStartTime = (trialStartIdx - 1) ./ srate;
trialStopTime  = (trialStopIdx - 1) ./ srate;

% Trial validity after artifact samples have been removed
validTrial = false(nTrials,1);
for tr = 1:nTrials
    validTrial(tr) = sum(trial_sc_vec == tr) >= length(window);
end

%% STFT across whole experiment
[STFTs, freqBins, t_ax] = stft(data, srate, ...
    'Window', window, ...
    'OverlapLength', noverlap, ...
    'FFTLength', nfft, ...
    'FrequencyRange', 'onesided');

% STFTs: [freq x time x channel]
nFreq = numel(freqBins);
nTime = numel(t_ax);

% Map STFT time-bin centers to samples so artifact-centered bins can be excluded.
stftSampleIdx = round(t_ax .* srate) + 1;
stftSampleIdx = max(1, min(numel(stim), stftSampleIdx));

%% Convert STFT coefficients to PSD
winNorm = srate * sum(window.^2);
STFT_psd = abs(STFTs).^2 ./ winNorm;

% One-sided correction
if rem(nfft,2) == 0
    if nFreq > 2
        STFT_psd(2:end-1,:,:) = 2 .* STFT_psd(2:end-1,:,:);
    end
else
    if nFreq > 1
        STFT_psd(2:end,:,:) = 2 .* STFT_psd(2:end,:,:);
    end
end

%% Mean PSD across whole experiment
artifactTimeMask = artifactMask(stftSampleIdx);
STFT_psd(:,artifactTimeMask,:) = NaN;
meanPSD_stft = squeeze(mean(STFT_psd, 2, 'omitnan')).';  % [channels x freq]
[meanPSD_welch, freq_bins] = pwelch(data(~artifactMask,:), window, noverlap, nfft, srate);

%% Baseline normalization
if isempty(baselinePSD)
    baseline = meanPSD_welch;
    baselineSource = "meanPSD";
else
    baseline = baselinePSD;
    if isequal(size(baseline), [nFreq nCh])
        baseline = baseline;
    elseif isequal(size(baseline), [nCh nFreq])
        baseline = baseline.';
    elseif ~isequal(size(baseline), [nCh nFreq])
        error('baselinePSD must be [channels x frequencies] or [frequencies x channels].');
    end
    baselineSource = "baselinePSD";
end
baseline(baseline == 0) = eps;

%% Initialize outputs
trialPSD = nan(nFreq, nCh,  nTrials);
trialNormPSD = nan(nFreq, nCh,  nTrials);
trialCPSD = cell(nTrials,1);
trialCoherence = cell(nTrials,1);

bandNames = fieldnames(bands);
bandPower = struct();

for b = 1:numel(bandNames)
    bandPower.(bandNames{b}) = nan(nCh, nTrials);
end

%% Trial-wise PSD, band power, CPSD, coherence
for tr = 1:nTrials
    if ~validTrial(tr)
        continue
    end

    % Select STFT time bins whose centers belong to clean samples in this trial
    trialTimeMask = trial_sc_vec(stftSampleIdx) == tr;

    if ~any(trialTimeMask)
        validTrial(tr) = false;
        continue
    end

    %% Trial PSD
    trial_data = data(trial_sc_vec == tr, :); % Artifact samples are excluded
    % Extract trial power spectrum (using pwelch)
    trialPxx =  pwelch(trial_data, window, noverlap, nfft, srate);

    % Alternatively from STFT PSD
    % trialPxx = squeeze(mean(STFT_psd(:,trialTimeMask,:), 2, 'omitnan')).'; 

    trialPSD(:,:,tr) = trialPxx;
    trialNormPSD(:,:,tr) = trialPxx ./ baseline;

    %% Band power
    for b = 1:numel(bandNames)
        currentBand = bandNames{b};
        fr = bands.(currentBand);
        fIdx = freqBins >= fr(1) & freqBins < fr(2);
        if any(fIdx)
            bandPower.(currentBand)(:,tr) = ...
                trapz(freqBins(fIdx), trialPxx(fIdx,:), 1);
        end
    end

    %% CPSD from STFT coefficients
    cpsdMatrix = complex(nan(nCh, nCh, nFreq));

    for ch1 = 1:nCh
        X1 = STFTs(:,trialTimeMask,ch1); % [freq x time]
        for ch2 = ch1:nCh
            X2 = STFTs(:,trialTimeMask,ch2);
            cpsd_f = mean(X1 .* conj(X2), 2, 'omitnan') ./ winNorm;
            % Same one-sided correction as PSD
            if rem(nfft,2) == 0
                if nFreq > 2
                    cpsd_f(2:end-1) = 2 .* cpsd_f(2:end-1);
                end
            else
                if nFreq > 1
                    cpsd_f(2:end) = 2 .* cpsd_f(2:end);
                end
            end
            cpsdMatrix(ch1,ch2,:) = cpsd_f;
            if ch1 ~= ch2
                cpsdMatrix(ch2,ch1,:) = conj(cpsd_f);
            end
        end
    end

    trialCPSD{tr} = cpsdMatrix;

    %% Coherence from CPSD
    trialCoherence{tr} = cpsdToCoherenceLocal(cpsdMatrix);
end

%% Package output
trialSpectFeatures = struct();

trialSpectFeatures.trialPSD = trialPSD;
trialSpectFeatures.meanPSD = meanPSD_welch;
trialSpectFeatures.meanPSD_stft = meanPSD_stft;
trialSpectFeatures.trialNormPSD = trialNormPSD;

trialSpectFeatures.TFmap = STFT_psd;
trialSpectFeatures.TFtime = t_ax;

trialSpectFeatures.validTrial = validTrial;
trialSpectFeatures.trialVec = trial_sc_vec; 
trialSpectFeatures.trialStimCode = tr_sc; 
trialSpectFeatures.trialStartIdx = trialStartIdx;
trialSpectFeatures.trialStopIdx = trialStopIdx;
trialSpectFeatures.trialStartTime = trialStartTime;
trialSpectFeatures.trialStopTime = trialStopTime;

trialSpectFeatures.freqBins = freqBins;
trialSpectFeatures.bandPower = bandPower;
trialSpectFeatures.trialCPSD = trialCPSD;
trialSpectFeatures.trialCoherence = trialCoherence;

trialSpectFeatures.info = struct();
trialSpectFeatures.info.srate = srate;
trialSpectFeatures.info.window = window;
trialSpectFeatures.info.noverlap = noverlap;
trialSpectFeatures.info.nfft = nfft;
trialSpectFeatures.info.bands = bands;
trialSpectFeatures.info.baselineSource = baselineSource;
trialSpectFeatures.info.minTrialSamples = length(window);
trialSpectFeatures.info.artifactCode = -1;

end

%% ------------------------------------------------------------------------
function validateBands(bands)

if ~isstruct(bands) || isempty(fieldnames(bands))
    error('bands must be a non-empty struct.');
end

bandNames = fieldnames(bands);

for b = 1:numel(bandNames)
    thisBand = bands.(bandNames{b});

    if ~isnumeric(thisBand) || numel(thisBand) ~= 2
        error('Band "%s" must be numeric [fLow fHigh].', bandNames{b});
    end

    if any(~isfinite(thisBand)) || any(thisBand < 0)
        error('Band "%s" must contain finite non-negative values.', bandNames{b});
    end

    if thisBand(1) >= thisBand(2)
        error('Band "%s" must satisfy fLow < fHigh.', bandNames{b});
    end

    if ~isvarname(bandNames{b})
        error('Band name "%s" is not a valid MATLAB field name.', bandNames{b});
    end
end
end

%% ------------------------------------------------------------------------
function cohMatrix = cpsdToCoherenceLocal(cpsdMatrix)

nCh = size(cpsdMatrix,1);
nFreq = size(cpsdMatrix,3);

cohMatrix = nan(nCh, nCh, nFreq);

autoPSD = zeros(nCh, nFreq);

for ch = 1:nCh
    autoPSD(ch,:) = real(squeeze(cpsdMatrix(ch,ch,:))).';
end

for ch1 = 1:nCh
    for ch2 = 1:nCh
        denom = autoPSD(ch1,:) .* autoPSD(ch2,:);
        cpsdVec = squeeze(cpsdMatrix(ch1,ch2,:)).';

        coh = nan(1,nFreq);
        valid = denom > 0;

        coh(valid) = abs(cpsdVec(valid)).^2 ./ denom(valid);
        coh(valid) = min(coh(valid), 1);

        cohMatrix(ch1,ch2,:) = reshape(coh,1,1,[]);
    end
end
end