function out = extractBSEPs_gtec(fpath, t_lims, keepChnls, rerefType)
% extractBSEPs_gtec
% Extract stimulation-evoked potentials from g.tec / BCI2000 recordings
%
% INPUTS
%   fpath      : full path to .dat file
%   t_lims     : [t_pre t_post] in seconds (default [-0.5 2])
%   keepChnls  : vector of channel indices to keep (default: all)
%
% OUTPUT
%   out        : struct containing stimulation responses and metadata

%% -------------------- Input handling --------------------
if nargin < 2 || isempty(t_lims)
    t_lims = [-0.5 2];
end

if nargin < 3
    keepChnls = [];
end

if nargin < 4
    rerefType = 'native';
end

%% -------------------- Load data --------------------
[signal, states, parameters] = load_bcidat(fpath, '-calibrated');
signal = double(signal);

nSamples = size(signal,1);
nCh      = size(signal,2);

if isempty(keepChnls)
    keepChnls = 1:nCh;
end

%% -------------------- File ID --------------------
[~, fileID] = fileparts(fpath);

%% -------------------- Sampling rate --------------------
srate = parameters.SamplingRate.NumericValue;

%% ---------------- Initiate output struct --------------------
out = struct();
out.FileID = fileID;

out.RecordingParameters = struct( ...
    'AmpGain', parameters.SourceChGain.Value{1}, ...
    'ReferenceChannel', parameters.RefChList.NumericValue, ...
    'FS', srate);

%% -------------------- Preprocessing --------------------
 % Hardcoded for the 1xdirectional and two 4x2 ECoG grids, future implementation will automatize
 % split into 3 leads (8 contacts each)
leads = {[1:8]; ...
             [9:12;13:16];...
             [17:20;21:24]};
cnstrct = 'combined'; spec = 'D-G-G-G';
bad_chnls = 25;

data = signal(:, keepChnls);
% Re-reference the data
switch rerefType
    case 'native'
        chNames = arrayfun(@(ch) sprintf('Ch%d', ch), keepChnls, 'UniformOutput', false);

    case 'DP' % bipolar within each lead
        % Hardcoded for the 4x2 ECoG grids
        % Future implementation of automatic re-refernce method will be here
        reref_data = [];
        chNames = {};
        [reref_data, chNames] = custom_reref(data,cnstrct,bad_chnls,'BP',spec,false);
        data = reref_data;
        clear reref_data
end
out.RecordingParameters.LeadOrganization = leads;
out.Preprocessing.ReRefType = rerefType;
out.Preprocessing.ChannelLabels = chNames;

% ~~~~~~~~~~~ Filter the data ~~~~~~~~~~~~~
disp('Filtering the data')
% High pass
% F_cutoff = 1; HP_order = 4;
% data = ieeg_highpass(data, srate);        % ieeG filter bank
% data = fastButterFilt(signal, FS, 'high', F_cutoff, HP_order); % Fast filter using CNEL mex multithread filter

%Notch
notchFreq = 60;
data = ieeg_notch(data, srate, 60);       % ieeG filter bank
% data = fastButterFilt(data, srate, 'notch', notchFreq); % Fast filter using CNEL mex multithread filter

%Save preprocessing paramters( NEED TO AUTOMATIZE, currrently hardcoded same as the filter implementations)
%from the ieeeg_highpass
% out.Preprocessing.HP.passFreq = 0.50;     % Hz
% out.Preprocessing.HP.stopFreq = F_cutoff;     % Hz
% out.Preprocessing.HP.passRipple= 3;        % dB
% out.Preprocessing.HP.stopAtten  = 30;       % dB

%NEED TO AUTOMATIZE, currrently hardcoded from the ieeeg_notch!
out.Preprocessing.Notch.FilterOrder = 4;
out.Preprocessing.Notch.FilterType  = 'butter';
out.Preprocessing.Notch.NotchFrequency = notchFreq;
out.Preprocessing.Notch.HalfPowerFrequencies  = [59 61]; %HARD CODED!!!!

%% -------------------- Stimulation parameters --------------------
stimParams = struct();
stimParams.Magnitude      = parameters.Magnitude;
stimParams.Modularity     = parameters.Modularity;
stimParams.Polarity       = parameters.Polarity;
stimParams.PhaseDuration  = parameters.PhaseDuration;

if isfield(parameters, 'NumberOfPulses')
    stimParams.NumberOfPulses     = parameters.NumberOfPulses.NumericValue;
    stimParams.FrequencyOfTrains  = parameters.FrequencyOfTrains.NumericValue;
    stimParams.NumberOfTrains     = parameters.NumberOfTrains.NumericValue;
end


%% -------------------- Parse stimulation pairs --------------------
stimConfig = parameters.StimConfig;
numPairs   = size(stimConfig.Value, 2);

stimPairs = zeros(numPairs, 2);
for k = 1:numPairs
    stimPairs(k,1) = str2double(stimConfig.Value{2,k}{1});
    stimPairs(k,2) = str2double(stimConfig.Value{3,k}{1});
end

%% -------------------- Identify stimulation epochs --------------------
StimulusCode = double(states.StimulusCode);
sc = zeros(size(StimulusCode));

for k = 1:numPairs
    a = stimConfig.Value{1,k};   % e.g., something like "StimulusCode==X&StimulusCode==Y"
    % --- Match the original code's string editing ---
    b = strfind(a,'&');
    if numel(b) >= 2
        a(b(1)) = ',';     % turn first & into a comma (for and(arg1,arg2))
        a(b(2)) = [];      % delete second &
        a(b(1)-1) = [];    % delete char just before first &, as in original script
    end

    % --- Evaluate the condition and get indices ---
    try
        eval(['idx=find(and(' a '));']);
    catch ME
        error('Failed to parse stimConfig expression for pair %d.\nExpression after edit:\n%s\nOriginal error:\n%s', ...
            k, a, ME.message);
    end

    % Reject very short stim segments (< 40 s)
    if numel(idx) / srate < 40
        sc(idx) = 0;
    else
        sc(idx) = k;
    end
end

%% -------------------- Select stimulation digital input --------------------
if max(states.DigitalInput1) == 0
    stimAll = double(states.DigitalInput2);
else
    stimAll = double(states.DigitalInput1);
end

%% -------------------- Time windowing --------------------
tIdx = round(t_lims(1)*srate) : round(t_lims(2)*srate);
tVec = tIdx / srate;

%% -------------------- Extract responses per stim pair --------------------
out.StimParameters = stimParams;
out.StimParameters.StimChannels.Source      = stimPairs(:,1);
out.StimParameters.StimChannels.Destination = stimPairs(:,2);
out.StimParameters.StimPeriodTimeLimits = t_lims;

out.StimResponse = struct();
out.StimResponse.Responses = cell(numPairs,1);
out.StimResponse.Time      = tVec;

for pair = 1:numPairs

    pairIdx = find(sc == pair);
    if isempty(pairIdx)
        continue
    end

    stim = stimAll(pairIdx);

    % Detect rising edges
    stimOnsets = find(diff(stim) == 1) + 1;

    % Remove onsets too close to end
    stimOnsets(stimOnsets + tIdx(end) > numel(pairIdx)) = [];

    if isempty(stimOnsets)
        continue
    end

    % Preallocate
    nStim = numel(stimOnsets);
    V = zeros(numel(tIdx), nStim, size(data,2));

    for k = 1:nStim
        V(:,k,:) = data(pairIdx(stimOnsets(k)) + tIdx, :);
    end

    out.StimResponse.Responses{pair} = V;
    out.StimResponse.StimPair{pair} = stimPairs(pair,:); 
    key = sprintf('StimPair_%d_%d', stimPairs(pair,1), stimPairs(pair,2));
    out.StimParameters.StimOnsets.(key) = stimOnsets;

end


end