clear all; close all; clc
%% Add paths if necesscary
bci2000path -AddAll
addpath('./functions/');
addpath('path/to/mnl_ieegBasics/functions'); % https://github.com/MultimodalNeuroimagingLab/mnl_ieegBasics

%% Select files and specify filepaths
data_folder = ['/path/to/data' filesep];
% %% Slect files manualy by specifyng paths
% fpth = '...';
% fnames = {'...',
%           '...'};

% -- Select runs with GUI  ----
[fnames, fpth] = uigetfile([data_folder '*.dat'],'Select functional .DAT FILE(s)','MultiSelect','on');
if isa(fnames,'char')   % prevents from getting error while loading data from single file
    fnames = {fnames};
end

% Sort Fnames so they are in ascending Session and Run 
tokens = regexp(fnames, 'S(\d+)R(\d+)', 'tokens'); % Extract session and run numbers
sesrun = cellfun(@(x) [str2double(x{1}{1}), str2double(x{1}{2})], ...
    tokens, 'UniformOutput', false); % Convert to numeric matrix [Session Run]
sesrun = vertcat(sesrun{:});
% Sort first by session, then by run
[~, idx] = sortrows(sesrun, [1 2]);
% Reorder filenames
fnames = fnames(idx);

clear sesrun tokens idx

%% Select baseline recording
[base_name, base_pth] = uigetfile([data_folder '*.dat'],'Select baseline .DAT FILE','MultiSelect','off');
if isa(fnames,'char')   % prevents from getting error while loading data from single file
    base_name = string(base_name);
    base_pth = string(base_pth);
end

%%  Specify Subject ID
subID = 'IMG';

%% Specify Constants and parameters
EMG_chans = struct( ...
        'hand_idx', [], ...
        'tongue_idx',[], ...
        'foot_idx',[]);
AMP_FACTORS = {'57.5dB','51.5dB','45.5dB','39.5dB'};    % set list of amplifaction factors for Cortec
interpMethod = 'linear';

%% Specify specral calculation parameters
% Specify frequency bands for statistical analysis
powerBands = struct('HFB', [70 110],...     % Specify broadband frequency range in Hz
    'LFB', [16 30]);      % Specify low-frequency power band in Hz

%% Specify Filtering parameters
% Select data filtration - if no filtration just comment out
BPf = [0.25 490];  % Specify Band Pass cut off frequency in Hz
BP_filt =true; % Specify whether band-pass filter should be used 
notch = true; % Specify whether notch filter should be used 
notchFreq = 60; 

%% Specify Spectral calculation parameters
SR  = 1000;  % !!! Hardcoded for specifyng spectral calclulation parameters !!!!
wndw_fxn = hann(2^nextpow2(0.5*SR), 'periodic');        % Hann
noverlap = floor(length(wndw_fxn) * 0.85);
nfft = 2*length(wndw_fxn);

%% Initialize variables
data = struct();

%% Load and extract data from baseline
[signal, states, parameters, total_samples ] ...      
        = LoadBCI2kData([base_pth base_name], interpMethod);  
signal = double(signal);
temp_name = strsplit(base_name,'.');
temp_name = matlab.lang.makeValidName(temp_name{1}); % Generate a valied struct fieldname
SR = parameters.SamplingRate.NumericValue;              % Get sampling frequency (rate)

transmitCh = parameters.TransmitChList.Value;
transmitCh = cellfun(@(x) str2double(regexp(x, '\d+', 'match', 'once')), transmitCh)'+1;

ref_sig = signal(:,transmitCh(1))-signal(:,transmitCh(2));

baselinePSD_raw = pwelch(signal,wndw_fxn,noverlap,nfft,SR);
[baselinePSD, freqBins] = pwelch(ref_sig,wndw_fxn,noverlap,nfft,SR);

data.baseline.FilePath = [base_pth temp_name];    % Store file name
data.baseline.FileName = temp_name; 
data.baseline.RawSignalPSD = baselinePSD_raw;
data.baseline.RerefSignalPSD = baselinePSD;  % Store rereferenced signal;
data.baseline.RerefSignalLabels = arrayfun(@(ch) sprintf('Ch%d', ch), transmitCh, 'UniformOutput', false);
data.baseline.SpectralInfo = struct('srate', SR, 'window', wndw_fxn, ...
    'noverlap', noverlap, 'nfft', nfft, 'freqBins', freqBins);
clear base_pth  base_name temp_name signal SR transmitCh
%% Load data & extract parameters from functional files
% Iterate through the file paths
for i = 1:length(fnames)
    disp(['Processing file ' fnames{i}])
    %% Load data
    [signal, states, parameters, total_samples ] ...      
        = LoadBCI2kData([fpth fnames{i}], interpMethod);  
    signal = double(signal);

    % Extract name
    fname_og =strsplit(fnames{i},'.');                         % Extract the filename without extension
    fname = matlab.lang.makeValidName(fname_og{1}); % Generate a valied struct fieldname

    % Extract parameters
    SR = parameters.SamplingRate.NumericValue;              % Get sampling frequency (rate)
    ref_ch = parameters.ReferenceCh.NumericValue;           % Get reference channel
    amp_f = AMP_FACTORS{parameters.AmplificationFactor.NumericValue + 1}; % Get the amplification factor

    % Additional parameters / variables
    s_fn = fieldnames(states);
    filterChain = parameters.SignalProcessingFilterChain.Value;
    t_ax = linspace(0,size(signal,1)/SR, size(signal,1));   % Creating time axis for plotting signal
    
    %% If SGNL states do not exist run BC2000 chain or if they are still in unit32
        sigSrcFilters = 'TransmissionFilter';    %strjoin(parameters.SignalSourceFilterChain.Value(hlp,1), '|');
        isSGNL = startsWith(s_fn, 'SGNL'); %   1) any field starting with 'SGNL'
        isNorm = ~cellfun(@isempty, regexp(s_fn, '^Normalizer(?:Offset|Gain)\d+$', 'once'));         %   2) NormalizerOffset%d or  NormalizerGain%d
        fieldsToConvert = s_fn(isSGNL | isNorm);
        if isempty(isSGNL) || isempty(isNorm)
            sigProcFilters = strjoin(parameters.SignalProcessingFilterChain.Value(1:end-1,1), '|'); % No control state filter
            filterChain = [sigSrcFilters '|' sigProcFilters];
            output = bci2000chain([fpth fnames{i}], filterChain);
            % Needs a bit of work to resample the output to the same size as
            % original signal
        end
    % Housekeeping
    clear isSGNL isNorm fieldsToConvert 


    %% Processs the signal
    % Transmited channels for online-pocessig (account for 0-based indexing)
    transmitCh = cellfun(@(s) str2double(regexp(s, '\d+$', 'match', 'once')), parameters.TransmitChList.Value)+1;
    
    % Save only the transmited channels
    raw_signal = signal(:,transmitCh);

    % Re-reference data
    [ref_sig, outputLabels, rerefInfo] = applyBCI2kSpatialFilter(signal(:,transmitCh), parameters, parameters.TransmitChList.Value); 

    if notch
        % ref_sig = ieeg_notch(ref_sig, SR, notchFreq, stopWidth); % Apply notch filter
        ref_sig = fastButterFilt(ref_sig, SR, 'notch', notchFreq, 4, true);
    end

    if BP_filt
        % signal = ieeg_butterpass(signal, BP_freqs, FS, true); % High-pass filter data using 3-ord Butter
        ref_sig = fastButterFilt(ref_sig, SR, 'bandpass', BPf, 3, true);
    end

    %% If EMG channel recordings are present save them
    if ~all(structfun(@isempty, EMG_chans)) % if at least one EMG channel is specified
        EMG = extractEMGs(signal, EMG_chans, SR, fname);
    end

    %% Catalogue trials
    bci_task = catalogue_BCI_trials(states, SR, 0.01, false);
    bci_task.Feedback = double(states.Feedback);            % For any other task
    bci_task.CursorPosY = double(states.CursorPosY);
    bci_task.TargetCode = double(states.TargetCode);
    bci_task.ResultCode = double(states.ResultCode);

    %%  Calculate spectograms and extract spectral feautres from each trial
    bci_control = bci_task.BCIcontrol;

    % Extract Spectral features
    spectralFeatures = trialSpectralFeatures(ref_sig, bci_control, SR, wndw_fxn, noverlap, nfft, powerBands, baselinePSD);
    
    % Nulling the notch filter dip in spectra
    % If not desired, just comment out the section 
    % ↓↓↓↓↓↓↓↓↓ from here ↓↓↓↓↓↓↓↓↓↓↓↓↓↓
    tps = spectralFeatures.TFmap;
    Pdb = pow2db(tps + eps);
    if notch
        help_idx = find((f<=notchFreq+1) & (f>=notchFreq-1));
        tps(help_idx,:,:) = nan;
        Pdb(help_idx,:,:) = nan;
    end
    %  ↑ ↑ ↑ ↑ ↑ to here  ↑ ↑ ↑ ↑ ↑ ↑ ↑ ↑
  
    clear help_idx
    
    %% Store the data
    % Store parameters
    data.subjectID = subID;
    data.(fname).FilePath = [fpth fnames{i}];    % Store file name
    data.(fname).FileName = fname_og; 
    data.(fname).SamplingRate = SR;         % Store Sampling Rate 
    data.(fname).Reference = ref_ch;        % Store reference channel
    data.(fname).AmpF = amp_f;              % Store amplification factor
    data.(fname).FilterChain = filterChain;   % Store filter chain
    if ~strcmp(parameters.TransmitChList.Value{1,1},'*') % Store transmited channels
        [~, data.(fname).ProccessedChanelsIdx] = ismember(parameters.TransmitChList.Value,...
            parameters.ChannelNames.Value);
    end
    data.(fname).FilterParameters = struct('Notch', notch, 'NotchFrequency', notchFreq, 'NotchStopWidth',1, 'NotchOrder',4,...
        'BandPass', BP_filt, 'BandPassFrequencies', BPf, 'BPorder', 3,...
        'FilterType', 'butter');  

    % Store signals
    data.(fname).RawSignal.Signal = raw_signal; % Store filtered signal Raw signal of transmited channels 
    data.(fname).RawSignal.ChannelLabels = arrayfun(@(ch) sprintf('Ch%d', ch), transmitCh, 'UniformOutput', false);
    data.(fname).RerefSignal.Signal = ref_sig;  % Store rereferenced signal;
    data.(fname).RerefSignal.SpectralFeatures = spectralFeatures; % Store extracted spectral features
    data.(fname).RerefSignal.Info = rerefInfo;

    % Store the info about perfomance
    data.(fname).BCITask = bci_task;

    % Store accelerometer/EMG data if they are present
    if exist('EMG','var')
        data.(fname).EMG = EMG;     
    end

    % If any BCI2000 outputs were saved into the states
    if any(startsWith(fieldnames(states), 'SGNL_')) 
        fldnames = fieldnames(states);
        matches = fldnames(startsWith(fldnames, 'SGNL_'));
        for fld = 1:numel(matches)
            data.(fname).BCI2000ProcessingOutputs.(matches{fld}) = typecast(states.(matches{fld}),'single');
        end
    end
    data.(fname).BCI2000ProcessingOutputs.classifierBins = parseBCI2000Classifier(parameters.Classifier);

    % Extract addiional parameters 
    if ismember('SpatialFilter',filterChain(:,1))
        data.(fname).SpatialFilter = parameters.SpatialFilterType;
    elseif ismember('LPFilter',filterChain(:,1))
        data.(fname).LPFilter = parameters.LPTimeConstant;
    end
    
end

%% Save 
outSave=questdlg(sprintf('Save MAT File and associated files in the same folder?\n %s',fpth), ...
    'Save Location for MAT File', ...
    'Yes','Choose location',['Don' '''' 't save' ],['Yes']);
switch outSave
    case 'Yes'
        saveloc = fpth;
        save_flag = true;
    case 'Choose location'
        saveloc = uigetdir('Choose a folder',fpth);
        save_flag = true;
    case ['Don' '''' 't save']
        save_flag = false;
end

if save_flag == true
    save([fpth subID '_BCI_data'], 'data');
elseif save_flag == false
    disp('File will NOT be saved')
else
    error('BCI2000 dat2mat CNEL save: Could not save file %s',filename);
end


