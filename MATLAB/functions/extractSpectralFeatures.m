function [spectralFeatures, STFTs, freqVec, t_ax, STFT_psd] = extractSpectralFeatures(signal, fs, window, overlap, nfft, bandsDefinition)
% extractSpectralFeatures
% Extracts spectral features from multichannel signal:
%   1) Band-power time courses
%   2) PSD per channel
%   3) CPSD matrix
%   4) Coherence matrix
%   5) Band-averaged coherence
%
% INPUTS
%   signal          : [samples x channels] or [channels x samples]
%   fs              : sampling frequency
%   window          : analysis window vector (default Hann, ~2 s, nextpow2)
%   overlap         : overlap in samples (default 50% of window)
%   nfft            : FFT length (default = window length)
%   bandsDefinition : optional struct with fields:
%                       bandsDefinition.bandName = [fLow fHigh]
%
% OUTPUT
%   spectralFeatures : struct with fields
%       .BandPower.(band)          [nTime x nCh]
%       .PSD                       [nFreq x nCh]
%       .CPSD                      [nCh x nCh x nFreq]
%       .Coherence                 [nCh x nCh x nFreq]
%       .BandCoherence.(band)      [nCh x nCh]
%       .info
% 

    % Defaults
    if nargin < 6 || isempty(bandsDefinition)
        bands = struct( ...
            'delta', [1 4], ...
            'theta', [4 8], ...
            'alpha', [8 15], ...
            'beta',  [15 30], ...
            'gamma', [30 100]);
    else
        % Validate custom bands
        if ~isstruct(bandsDefinition) || isempty(fieldnames(bandsDefinition))
            error('bandsDefinition must be a non-empty struct with fields like bandsDefinition.alpha = [8 15].');
        end
        bandNames = fieldnames(bandsDefinition);
        for b = 1:numel(bandNames)
            thisBand = bandsDefinition.(bandNames{b});
            if ~isnumeric(thisBand) || numel(thisBand) ~= 2
                error('Band "%s" must be a numeric 1x2 vector [fLow fHigh].', bandNames{b});
            end
            if any(~isfinite(thisBand)) || any(thisBand < 0)
                error('Band "%s" must contain finite, non-negative frequencies.', bandNames{b});
            end
            if thisBand(1) >= thisBand(2)
                error('Band "%s" must satisfy fLow < fHigh.', bandNames{b});
            end
            if ~isvarname(bandNames{b})
                error('Band name "%s" is not a valid MATLAB field name.', bandNames{b});
            end
        end
        bands = bandsDefinition;
    end

    if nargin < 5 || isempty(nfft)
        nfft = length(window);
    end
    if nargin < 4 || isempty(overlap)
        overlap = floor(length(window)/2);
    end
    if nargin < 3 || isempty(window)
        if nargin < 2 || isempty(fs)
            fs = 1000;
        end
        window = hann(2^(nextpow2(fs*2)), 'periodic');
    end
    if nargin < 2 || isempty(fs)
        error('Specify sampling frequency')
    end

    % Ensure [timeSamples x channels]
    if size(signal,2) > size(signal,1)
        signal = signal';
    end

    % Validate
    if nfft < length(window)
        error('nfft must be greater than or equal to the window length.');
    end
    if overlap >= length(window)
        error('overlap must be smaller than the window length.');
    end

    spectralFeatures = struct();


    % STFT: output dims are [freq x time x channel] with time across columns by default
    [STFTs, freqVec, t_ax] = stft(signal, fs, ...
        'Window', window, ...
        'OverlapLength', overlap, ...
        'FFTLength', nfft, ...
        'FrequencyRange', 'onesided');

    nFreq = numel(freqVec);
    nCh   = size(signal, 2);
    nTime = size(STFTs, 2);

    % Normalize to PSD-like units
    winNorm = fs * sum(window.^2);
    STFT_psd = abs(STFTs).^2 / winNorm;   % [freq x time x ch]

    % One-sided correction: double interior bins
    if nfft > 1
        if rem(nfft, 2) == 0
            % even nfft: keep DC and Nyquist unchanged
            if nFreq > 2
                STFT_psd(2:end-1,:,:) = 2 * STFT_psd(2:end-1,:,:);
            end
        else
            % odd nfft: keep DC unchanged, no exact Nyquist bin
            if nFreq > 1
                STFT_psd(2:end,:,:) = 2 * STFT_psd(2:end,:,:);
            end
        end
    end

    % ---------------------------------------------------------------------
    % 1) Band power time development 
    % ---------------------------------------------------------------------
    % Note: On May 5, 2026 (~10PM), band-power calculation was changed from
    % mean PSD within the specified band (units: µV^2/Hz) to integrated
    % band power over frequency (units: µV^2).

    bandNames = fieldnames(bands);
    for b = 1:numel(bandNames)
        currentBand = bandNames{b};
        fr = bands.(currentBand);
        fIdx = freqVec >= fr(1) & freqVec < fr(2);

        if any(fIdx)
            % Integrated band power from PSD -> [nTime x nCh]
            powerInBand = squeeze(trapz(freqVec(fIdx), STFT_psd(fIdx,:,:), 1));
            if nCh == 1
                powerInBand = powerInBand(:); % [nTime x 1]
            end
        else
            powerInBand = nan(nTime, nCh);
        end

        spectralFeatures.BandPower.(currentBand) = powerInBand;
    end

    % ---------------------------------------------------------------------
    % 2) CPSD matrix
    % ---------------------------------------------------------------------
    cpsdMatrix = complex(nan(nCh, nCh, nFreq));

    for ch1 = 1:nCh
        X1 = STFTs(:,:,ch1); % [freq x time]
        for ch2 = ch1:nCh
            X2 = STFTs(:,:,ch2); % [freq x time]
            cpsd_f = mean(X1 .* conj(X2) / winNorm, 2); % [freq x 1]

            % Apply same one-sided correction to CPSD
            if nfft > 1
                if rem(nfft, 2) == 0
                    if nFreq > 2
                        cpsd_f(2:end-1) = 2 * cpsd_f(2:end-1);
                    end
                else
                    if nFreq > 1
                        cpsd_f(2:end) = 2 * cpsd_f(2:end);
                    end
                end
            end

            cpsdMatrix(ch1, ch2, :) = cpsd_f;
            if ch1 ~= ch2
                cpsdMatrix(ch2, ch1, :) = conj(cpsd_f);
            end
        end
    end

    % ---------------------------------------------------------------------
    % 3) Auto-PSD from CPSD diagonal
    % ---------------------------------------------------------------------
    autoPSD = zeros(nFreq, nCh);
    for ch = 1:nCh
        autoPSD(:, ch) = real(squeeze(cpsdMatrix(ch, ch, :)));
    end

    % ---------------------------------------------------------------------
    % 4) Coherence
    % ---------------------------------------------------------------------
    cohMatrix = nan(nCh, nCh, nFreq);
    autoPSD_row = autoPSD.'; % [nCh x nFreq]

    for ch1 = 1:nCh
        for ch2 = 1:nCh
            denom = autoPSD_row(ch1,:) .* autoPSD_row(ch2,:);
            numer = abs(squeeze(cpsdMatrix(ch1, ch2, :))).'.^2;

            coh = nan(1, nFreq);
            valid = denom > 0;
            coh(valid) = numer(valid) ./ denom(valid);
            coh(valid) = min(coh(valid), 1); % numerical safety

            cohMatrix(ch1, ch2, :) = coh;
        end
    end

    % ---------------------------------------------------------------------
    % 5) Band-averaged coherence
    % ---------------------------------------------------------------------
    for b = 1:numel(bandNames)
        currentBand = bandNames{b};
        fr = bands.(currentBand);
        fIdx = freqVec >= fr(1) & freqVec < fr(2);

        if any(fIdx)
            spectralFeatures.BandCoherence.(currentBand) = mean(cohMatrix(:,:,fIdx), 3, 'omitnan');
        else
            spectralFeatures.BandCoherence.(currentBand) = ...
                nan(nCh, nCh);
        end
    end

    % Main outputs
    spectralFeatures.PSD = autoPSD;                 % [nFreq x nCh]
    spectralFeatures.CPSD = cpsdMatrix;             % [nCh x nCh x nFreq]
    spectralFeatures.Coherence = cohMatrix;         % [nCh x nCh x nFreq]

    % Save metadata
    spectralFeatures.info.nfft = nfft;
    spectralFeatures.info.window = window;
    spectralFeatures.info.overlap = overlap;
    spectralFeatures.info.SamplingFrequency = fs;
    spectralFeatures.info.Frequencies = freqVec;
    spectralFeatures.info.TimePoints = t_ax;
    spectralFeatures.info.Bands = bands;
end