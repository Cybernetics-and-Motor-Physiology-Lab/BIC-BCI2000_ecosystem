% FXNL_PIPELINE_01_EXTRACTDATA
% Functional-data extraction and preprocessing pipeline for BCI2000 recordings.
%
% Author: Frederik Lampert
% Institution: Mayo Clinic
% Year: 2026
%
% This script loads BCI2000 .dat recordings, applies optional signal filtering,
% extracts stimulation and task-related states, identifies bad channels, and
% stores the processed data in a MATLAB structure.
%
% Data availability:
%   Source recordings are available in BCI2000 .dat format through OpenNeuro:
%   doi:10.18112/openneuro.ds004624.v3.0.0
%
% Requirements:
%   - MATLAB R2023b
%   - BCI2000 MATLAB MEX tools (including load_bcidat)
%   - mnl_ieegBasics:
%       https://github.com/MultimodalNeuroimagingLab/mnl_ieegBasics
%   - Project-specific helper functions used below, including:
%       interpolate_lost_samples, runlength, identifyBadChannels,
%       segment_signal, ieeg_butterpass, ieeg_butterhighpass,
%       ieeg_butterlowpass, and ieeg_notch
%
% Before running:
%   1. Replace the placeholder paths below with local paths.
%   2. Select the appropriate subject in the "Specify constants" section.
%   3. Review filtering parameters and subject-specific bad-channel settings.
%
clear all; close all; clc
%% Configure paths
% Replace these placeholders with local paths before running.
projectPath = '/path/to/project/';
ieegBasicsPath = '/path/to/mnl_ieegBasics/functions/';
colormapPath = '/path/to/colormaps/';
data_folder = ['/path/to/data/' filesep];

addpath(projectPath);
addpath(ieegBasicsPath);
addpath(genpath(colormapPath));

%% Select files and specify file paths

% -- Select Functional Data with GUI  ----
[fnames, fpth] = uigetfile([data_folder '*.dat'], ...
    'Select functional .DAT FILE(s)', 'MultiSelect', 'on');

if isequal(fnames, 0)
    error('No BCI2000 .dat files were selected.');
elseif ischar(fnames)
    fnames = {fnames};  % Normalize single-file selection to a cell array
end


%% Specify constants 
sub = 'c01'; 

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
    case 'c05'
        cnstrct = 'combined'; spec = 'G-G-G-D';
        knownBadChnls = [];
end


AMP_FACTORS = {'57.5dB','51.5dB','45.5dB','39.5dB'};    % set list of amplification factors for Cortec
bCh_thr = 3; % Bad channel detection threshold (in #std)

% Select data filtering - if no filtering just comment out
HP = 0.5; % Specify High Pass cut off frequency in Hz
% LP = 350; % Specify Low Pass cut off frequency in Hz
notch = false; % Specify whether notch filter should be used

%% Initialize variables
data = struct();

%% Load data & extract parameters
% Iterate through the file paths
for i = 1:length(fnames)
    % Load data
    [signal, states, parameters, total_samples ] ...      
        = load_bcidat([fpth fnames{i}], '-calibrated');         % load using MEX package   
    signal = interpolate_lost_samples(signal,states,'linear');  % linear interpolation of lost samples
    
    % Extract filename
    fname =strsplit(fnames{i},'.');                         % Extract the filename without extension
    temp_name = replace(fname{1},'-','_'); % Generate a valid structure field name

    % Extract parameters
    SR = parameters.SamplingRate.NumericValue;              % Get sampling frequency (rate)
    ref_ch = parameters.ReferenceCh.NumericValue;           % Get reference channel
    amp_f = AMP_FACTORS{parameters.AmplificationFactor.NumericValue + 1};
    t_ax = linspace(0,size(signal,1)/SR, size(signal,1));   % Creating time axis for plotting signal
    fldnames = fieldnames(states);
    filterChain = parameters.SignalProcessingFilterChain.Value;

    % Filter the data
    if exist('HP') && exist('LP') 
        filt_sig = ieeg_butterpass(signal, [HP LP], SR);
        data.(temp_name).Preprocessing.HP = HP;     % Store HP parameters
        data.(temp_name).Preprocessing.LP = LP;     % Store LP parameters
    elseif exist('HP')
        filt_sig = ieeg_butterhighpass(signal,HP,SR);
        data.(temp_name).Preprocessing.HP = HP;     % Store HP parameters
    elseif exist('LP')
        filt_sig = ieeg_butterlowpass(signal,LP,SR);
        data.(temp_name).Preprocessing.LP = LP;     % Store LP parameters
    else 
        filt_sig = signal;
    end
    
    if notch        % If notch  filter is required
        filt_sig = ieeg_notch(filt_sig, SR, 60, [59 61]);
        data.(temp_name).Preprocessing.Notch = notch;  % Store if notch was used parameters
    end

    % % Extract Stimulus Code
    sc = double(states.StimulusCode);            % For any other task
    % sc = double(states.PhaseInSequence==2);        % For visual stimulation (flickering) task
    sc_dur = runlength(sc, total_samples);
    sc_onsets = [1; cumsum(sc_dur(1:end-1)) + 1];

    % Store the data
    data.(temp_name).FilePath = [fpth fnames{i}];    % Store file name
    data.(temp_name).FileName = fname{1};            % File name without extension
    data.(temp_name).Signal = filt_sig;    % Store filtered signal
    data.(temp_name).StimulusCode = sc;    % Store stimulus Code
    data.(temp_name).SequenceOnset = sc_onsets;
    data.(temp_name).SequenceDuration = sc_dur;

    data.(temp_name).BadChannels = identifyBadChannels(filt_sig, cnstrct, bCh_thr); % Suggested bad channels
    data.(temp_name).SamplingRate = SR;         % Store Sampling Rate 
    data.(temp_name).Reference = ref_ch;        % Store reference channel
    data.(temp_name).AmpF = amp_f;              % Store amplification factor
    data.(temp_name).FilterChain = filterChain;   % Store filter chain
    if ~strcmp(parameters.TransmitChList.Value{1,1},'*') % Store transmitted channels
        [~, data.(temp_name).ProcessedChannelsIdx] = ismember(parameters.TransmitChList.Value,...
            parameters.ChannelNames.Value);
    end
    

    % Extract the stimulation annotation
    if parameters.EnableStimulation.NumericValue == 1
        data.(temp_name).ImplantStimulation = single(states.ImplantStimulation);
        data.(temp_name).RequestedStimulation = single(states.RequestedStimulation);
    end

    % Extract additional states
    % Accelerometer and audio-envelope states are handled independently:
    %   - if both exist, both are extracted and passed to segment_signal;
    %   - if only one exists, only the available modality is extracted;
    %   - if neither exists, the original StimulusCode remains the task annotation.
    hasAcc = all(isfield(states, {'MTw1X','MTw1Y','MTw1Z'}));
    hasAudio = isfield(states, 'AudioInEnvelope1');

    if hasAcc || hasAudio
        annotationData = [];
        annotationLabels = {};

        if hasAcc
            accX = (double(states.MTw1X) - 2^31) ./ parameters.DataPrescaler.NumericValue;
            accY = (double(states.MTw1Y) - 2^31) ./ parameters.DataPrescaler.NumericValue;
            accZ = (double(states.MTw1Z) - 2^31) ./ parameters.DataPrescaler.NumericValue;
            accData = [accX accY accZ];

            data.(temp_name).AccData = accData;
            annotationData = [annotationData accData]; %#ok<AGROW>
            annotationLabels = [annotationLabels {'AccX','AccY','AccZ'}]; %#ok<AGROW>
        end

        if hasAudio
            audio = double(states.AudioInEnvelope1) ./ 100;

            data.(temp_name).AudioEnvelope = audio;
            annotationData = [annotationData audio]; %#ok<AGROW>
            annotationLabels = [annotationLabels {'AudioEnvelope'}]; %#ok<AGROW>
        end

        % Store the complete auxiliary annotation matrix and its column labels.
        data.(temp_name).AnnotatedData = annotationData;
        data.(temp_name).AnnotatedDataLabels = annotationLabels;

        % Launch the manual/interactive segmentation tool using whichever
        % auxiliary modalities are available.
        [marked_times, marked_ids, modalities] = ...
            segment_signal(annotationData, sc, SR, fname{1});

        % Generate a binary vector indicating intervals marked as task periods.
        taskStatus = zeros(size(signal,1),1);
        marked_idxs = round(marked_times .* SR);

        % Keep indices inside the recording boundaries.
        marked_idxs = max(1, min(size(signal,1), marked_idxs));

        for j = 1:2:(numel(marked_idxs)-1)
            taskStatus(marked_idxs(j):marked_idxs(j+1)) = 1;
        end

        data.(temp_name).Onsets = [marked_ids' marked_times'];
        data.(temp_name).Events = modalities;
        data.(temp_name).isPerformingTask = taskStatus;

    else
        % No accelerometer or audio-envelope states are available.
        % Keep the BCI2000 StimulusCode as the task annotation rather than
        % attempting auxiliary-state segmentation.
        data.(temp_name).isPerformingTask = sc;
    end

    % Extract other optional BCI2000 processing states independently of the
    % accelerometer/audio-envelope annotation above.
    if any(startsWith(fldnames, 'SGNL_'))
        matches = fldnames(startsWith(fldnames, 'SGNL_'));
        for fld = 1:numel(matches)
            data.(temp_name).BCI2000ProcessingOutputs.(matches{fld}) = ...
                typecast(states.(matches{fld}), 'single');
        end
    end

    if all(isfield(states, {'ControlState','RefractoryPeriod'}))
        data.(temp_name).ControlState.Vec = single(states.ControlState);
        data.(temp_name).ControlState.Expression = parameters.ControlExpression.Value;
        data.(temp_name).RefractoryPeriod.Vec = single(states.RefractoryPeriod);
        data.(temp_name).RefractoryPeriod.Value = parameters.RefractoryPeriod.Value;
    end
    % Extract additional parameters 
    if ismember('SpatialFilter', filterChain(:,1))
        data.(temp_name).SpatialFilter = parameters.SpatialFilterType;
    end
    if ismember('LPFilter', filterChain(:,1))
        data.(temp_name).LPFilter = parameters.LPTimeConstant;
    end
    
end
%% Save
outSave = questdlg( ...
    sprintf('Save MAT file and associated files in the same folder? %s', fpth), ...
    'Save Location for MAT File', ...
    'Yes', 'Choose location', 'Don''t save', 'Yes');

switch outSave
    case 'Yes'
        saveloc = fpth;
        save_flag = true;

    case 'Choose location'
        saveloc = uigetdir(fpth, 'Choose output folder');
        save_flag = ~isequal(saveloc, 0);

    otherwise
        save_flag = false;
end

if save_flag
    save(fullfile(saveloc, 'fxnl_data.mat'), 'data', '-v7.3');
    fprintf('Saved functional data to: %s', fullfile(saveloc, 'fxnl_data.mat'));
else
    disp('File will NOT be saved.')
end
