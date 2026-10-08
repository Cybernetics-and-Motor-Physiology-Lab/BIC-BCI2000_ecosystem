clear all; close all; clc
%% Add paths if necesscary
addpath(genpath('./functions'));
addpath(genpath('./colormaps'));
bci2000path -AddAll

% Dependencies:
%    - BCI2000 mx fucntions
%    - LoadBCI2kData  - for loading data 
%    - extractImpedanceValues
%    - identifyBadChannels
%    - runlength
%    - extractSpectralFeatures
%    - struct2single
%    - spikeDetector
%    - progressbar - not strictly necessary
%    - parseBCI2000_parameters_values

%% Specify file paths to raw files

baseline = 'path/to/baselineRec.dat';

stim_files = {
    '/path/to/StimData1.dat';...
    '/path/to/StimData2.dat';...
    '/path/to/StimData3.dat';...
    };


%% Specify Constants 
sub = 'c04'; 

switch sub
    case 'c01'
        cnstrct = 'surface'; spec = 'vertical';  % Custom difference pair for cortec Grids
        knownBadChnls= [13,22];
    case 'c02'
        cnstrct = 'depth'; spec = 'D-S-D-S';
        knownBadChnls= [1:8];
    case 'c03'
        cnstrct = 'surface'; spec = 'horizontal';  % Custom difference pair for cortec Grids
        knownBadChnls= [11,21];
    case 'c04'
        cnstrct = 'combined'; spec = 'G-G-D-S';
        knownBadChnls = 1:8;
        channelLeadMap = dictionary( ...
        'displaced', struct('indices', 1:8, 'labels',  "Ch" + (1:8)), ...
        'occipital',struct('indices', 9:16, 'labels',  "Ch" + (9:16)), ...
        'thalamus',struct('indices', 17:24, 'labels',  "Ch" + (17:24)),...
        'lobus-piriformis', struct('indices', 25:32, 'labels',  "Ch" + (25:32)) );
    case 'c05'
        cnstrct = 'combined'; spec = 'G-G-G-D';
        knownBadChnls = [];
end

AMP_FACTORS = {'57.5dB','51.5dB','45.5dB','39.5dB'};    % set list of amplifaction factors for Cortec

%% Specify data processing parameters 

% Specify initial data curation parameters
    bCh_thr = 3; % Bad channel detection threshold (in #std)
    interpMethod = 'linear'; % Specify interpolation method to replace lost samples

% Select spatial re-referencing
    % reref = 'None';
    % reref = 'CAR';         % Common average
    % reref = 'BP'; 
    % reref = 'leadSpecific'; RR_chnls = [1 9 17 25]
    reref = 'leadCAR';   % Custom Common Average for cortec Grids

% Select data filtration parameters
    BP_filt = true;
    BP_freqs = [0.5 499]; % Specify High Pass cut off frequency in Hz
    notch = false; % Specify whether notch filter should be used

% Specify spike detecion parameters
    SpikeDetectionParams = struct('LF',10, ...
        'HF',85, ...
        'th_mult', 3.65, ...
        'Amp_th_min',85);

% Specify STFT parameters
    wndw_length = 2^nextpow2(1000);
    spectralParameters = struct('window', wndw_length, ...
                'overlapCoeff',0.8, ...
                'windowFunction','hann', ...
                'NFFT', wndw_length, ...
                'bands', struct('Delta', [1 4], ...
                'Theta', [4 8], ...
                'Alpha', [8 15], ...
                'Beta',  [15 30], ...
                'Gamma', [30 100]));
    clear wndw_length

%% Initialize variables
data = struct();

%% Load and extract data 
files = [baseline; stim_files]; % Concatenate

for f=1:numel(files)
    [filepath,fname_og,ext] = fileparts(files(f));
    fname = matlab.lang.makeValidName(fname_og{1});
    fprintf('Processing file %s. (%d out of %d)\n', fname, f,  numel(files))

    [signal, states, parameters, total_samples ] ...      
            = LoadBCI2kData(files{f}, interpMethod);  

    % Load impeance values if they exist
    if isfile(replace(files(f),'.dat', '_info1.txt'))
        impedanceValues = extractImpedanceValues(replace(files{f},'.dat', '_info1.txt'));
    else
        impedanceValues = [];
    end

    FS = parameters.SamplingRate.NumericValue;              % Get sampling frequency (rate)

    %% Extract stim and rest periods 
    if parameters.EnableStimulation.NumericValue == 1
        sc_dur = [];
        sc = [];
        sc_onsets = [];
        [sc_dur, sc] = runlength(states.ImplantStimulation, total_samples);
        sc_onsets = [1; cumsum(sc_dur(1:end-1)) + 1];
        % Rest periods
        rest_idx = [sc_onsets(sc==0), sc_onsets(sc==0)+sc_dur(sc==0)-1];
        rest_indices = arrayfun(@(s, e) s:e, rest_idx(:,1), rest_idx(:,2), 'UniformOutput', false);
        % Stim periods
        stim_idx = [sc_onsets(sc==1), sc_onsets(sc==1)+sc_dur(sc==1)-1];
        stim_indices = arrayfun(@(s, e) s:e, stim_idx(:,1), stim_idx(:,2), 'UniformOutput', false);
    else
        rest_idx = [1 total_samples];
        rest_indices = {1:total_samples};
    end

    %% Identify bad channels
    [bad_chnls, badChannelInfo] = identifyBadChannels(signal([rest_indices{:}], :), cnstrct, FS,bCh_thr, knownBadChnls, impedanceValues, 0);
    g_chnls = 1:32; g_chnls(bad_chnls) = [];

    %% Spatial re-referncing
    switch reref
        case 'Native'
            proc_signal = signal;
            ch_names = arrayfun(@(ch) sprintf('Ch %d', ch), 1:size(signal,2), 'UniformOutput', false);
            ch_names(bad_chnls) = arrayfun(@(n) sprintf('BAD-Ch %d', n), bad_chnls, 'UniformOutput', false);
        case 'BP'
            [proc_signal, ch_names] = custom_reref(signal, cnstrct, bad_chnls, reref, spec, true);
        case 'leadSpecific'
            [proc_signal, ch_names] = custom_reref(signal, cnstrct, bad_chnls, reref, RR_chnls, true);
        otherwise
            [proc_signal, ch_names] = custom_reref(signal, cnstrct, bad_chnls, reref, [], true);
    end


    %% Filtering 
    if notch
        % signal = ieeg_notch(signal, FS, 60, stopWidth); % Apply notch filter
        proc_signal = fastButterFilt(proc_signal, FS, 'notch', 60, 4, true);
    end

    if BP_filt
        % signal = ieeg_butterpass(signal, BP_freqs, FS, true); % High-pass filter data using 3-ord Butter
        proc_signal = fastButterFilt(proc_signal, FS, 'bandpass', BP_freqs, 3, true);
    end

    %% Temporal features
    stats = struct();
    stats.Average = mean(proc_signal, 1);
    stats.Median = median(proc_signal, 1);
    stats.STD = std(proc_signal, 1);
    stats.Variance = var(proc_signal, 1);
    stats.RMS = rms(proc_signal, 1);
    stats.Skewness = skewness(proc_signal, 1);
    stats.Kurtosis = kurtosis(proc_signal, 1);
    stats.ZeroCrossRate = zerocrossrate(proc_signal);
    stats.AutoCorr = corr(proc_signal);
    stats.AutoCov = cov(proc_signal);

%% Initialize Spectral parameters
    switch lower(spectralParameters.windowFunction)
        case 'hann'
            wndw_fxn = hann(spectralParameters.window, 'periodic');
        case 'hamming'
            wndw_fxn = hamming(spectralParameters.window, 'periodic');
        otherwise
            error('Unsupported window function: %s', spectralParameters.windowFunction);
    end
    noverlap = floor(length(wndw_fxn) * spectralParameters.overlapCoeff);
    nfft = spectralParameters.NFFT;
    bands = spectralParameters.bands;

%% Initialize structs
    spectral_features = struct();
    temporal_features = struct();
    SpikeDetect = struct();
    rest_sig = [];

%% Feature extraction
    fprintf('Extracting features from %s. \n', fname)
    for s = 1:size(rest_idx,1)
        % --------------- Extract rest periods ---------
        rest_sig = proc_signal(rest_idx(s,1):rest_idx(s,2),:);
        trialName = sprintf('Trial%d', s);
    
        % --------------- Temporal features ----
        if size(rest_sig,1) < 2
            temporal_features.(trialName).Covariance = [];
            temporal_features.(trialName).CrossCorrelation = [];
            temporal_features.(trialName).MutualInformation = [];
        else
            temporal_features.(trialName).Covariance = single(cov(rest_sig, 'omitrows'));
            temporal_features.(trialName).CrossCorrelation = single(corrcoef(rest_sig, 'Rows','pairwise'));
            temporal_features.(trialName).MutualInformation = []; % This still needs to be implemented
        end
    
        % --------------- Spectral features -----
        if size(rest_sig,1) < length(wndw_fxn)
            spectral_features.(trialName).PSD = [];
            spectral_features.(trialName).Coherence = [];
        else
            % Caluclate and extract spectral features
            [spectralFeatures, ~, freqVec, t_ax, STFT_psd] = extractSpectralFeatures(rest_sig, FS, wndw_fxn, noverlap, nfft, bands);
            % Convert to singles and save
            spectral_features.(trialName) = struct2single(spectralFeatures);
        end
    
        % -------------------- Detect spikes ------------
        % Upsample signal to 2k Hz
        q = 2;
        rest_sig = resample(rest_sig,q,1,'Dimension',1);
        %Create montage
        montage = struct('SampleRate', q*FS, 'ChannelNames', {ch_names});
        %Detect spikes
        % For spike detection and elimination functions please refer to
        % Ayyoubi, A. H., Besheli, B. F., Swamy, C. P., Okkabaz, J. L., Miller, K. J., Worrell, G. A., & Ince, N. F. (2025).
        % Spurious Spike Elimination using Sparse Signal Processing Improves Seizure Onset Zone Delineation in Brief Intraoperative iEEG Recordings.
        % The ... Midwest Symposium on Circuits and Systems conference proceedings : MWSCAS. Midwest Symposium on Circuits and Systems, 2025, 666–670.
        % https://doi.org/10.1109/mwscas53549.2025.11244470
        % Filter spikes
        IIS_Det = Spike_detection_function(rest_sig, montage, [], 512, ...
            [SpikeDetectionParams.LF, SpikeDetectionParams.HF], ...
            SpikeDetectionParams.Amp_th_min, SpikeDetectionParams.th_mult);
        Spike_RFOMP = pSpike_elimination_RFOMP(IIS_Det.spike,montage.SampleRate);
        denoisedSpikesIdx = find(Spike_RFOMP.feature.Pred==1);  % [samples x nSpikes]
    
        SpikeDetect = IIS_Det.spike;
        SpikeDetect.Montage = montage; % !!!! changed from channel names
        SpikeDetect.RFOMP = Spike_RFOMP; % Save the output of spike filtering
        SpikeDetect.UpsamplingFactor = q;

        progressbar(s, size(rest_idx,1))  % Console progress, if causng problem just comment out
    end

    %% Pack everything into struture
    % General info
    data.(fname).folder = filepath{1};
    data.(fname).name = [fname_og{1} ext{1}];
    data.(fname).SamplingRate = FS;
    data.(fname).NativeRef = parameters.ReferenceCh.NumericValue;
    data.(fname).Duration = string(seconds(total_samples/FS),"hh:mm:ss.SSS");
    data.(fname).PacketLoss = calculatePacketLoss(states);
    data.(fname).AmplFactor = AMP_FACTORS{parameters.AmplificationFactor.NumericValue+1};
    data.(fname).Ground = parameters.UseGround.NumericValue;
    data.(fname).BadChannels.BadChannels = bad_chnls;
    data.(fname).BadChannels.Info = badChannelInfo;
    data.(fname).Impedance = impedanceValues;

    % Save proccessed signal
    data.(fname).Signal = proc_signal;

    % Preprocessing Info
    data.(fname).Preprocessing.LostSampleInterpolation = interpMethod;
    data.(fname).Preprocessing.Filtering.BandPass = BP_freqs;
    data.(fname).Preprocessing.Filtering.Notch = notch;
    data.(fname).Preprocessing.Rerefencing.Type = reref;
    data.(fname).Preprocessing.Rerefencing.Construct = cnstrct;
    data.(fname).Preprocessing.Rerefencing.ConstructSpecification = spec;
    data.(fname).ChannelInfo.ChannelNames = ch_names;
    data.(fname).ChannelInfo.LeadMap = channelLeadMap;

    % Stim info
    if parameters.EnableStimulation.NumericValue == 1
        data.(fname).StimInfo.StimPulses = parseBCI2000_parameters_values(parameters, 'StimulationPulses');
        data.(fname).StimInfo.StimTriggers = parseBCI2000_parameters_values(parameters, 'StimulationTriggers');
        data.(fname).StimInfo.StimSequence.Sequence = sc;
        data.(fname).StimInfo.StimSequence.SequenceOnset = sc_onsets;
        data.(fname).StimInfo.StimSequence.SequenceDuration = sc_dur;
    end

    % General stats
    data.(fname).GeneralStats = stats;
    % Temporal features
    data.(fname).TemporalFeatures = temporal_features;
    % Spectral Features
    data.(fname).SpectralFeatures = spectral_features;
    data.(fname).SpectralFeatures.TFmap.PSDs = STFT_psd;
    data.(fname).SpectralFeatures.TFmap.Time = t_ax;
    data.(fname).SpectralFeatures.TFmap.FreqBins = freqVec;
    % Spikes
    data.(fname).SpikeDetection = SpikeDetect;

    fprintf('Finished processing %s file. \n', fname)
end

%% Save 
[savefile, saveloc] = uiputfile('stim_data.mat');
save([saveloc savefile], 'data', '-v7.3');

