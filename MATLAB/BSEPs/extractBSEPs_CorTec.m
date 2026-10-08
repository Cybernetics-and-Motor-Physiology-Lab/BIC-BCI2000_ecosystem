% Script to load, extract, pre-process and reshape data  from CorTec
% measurements for CRP analysis

% Dependencies:
%    - BCI2000 mx fucntions
%    - LoadBCI2kData  - for loading data
%    - identifyBadChannels
%    - runlength
%    - parseBCI2000_parameters_values

%%
clear all; close all; clc
%% ~~~~~~~~~~~~~~~~~~~~ ADD PATHS ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
addpath('path/to/BCI2000/mex');
addpath('path/to/crp_scripts'); % https://github.com/kaijmiller/crp_scripts
addpath('path/to/mnl_ieegBasics/functions'); % https://github.com/MultimodalNeuroimagingLab/mnl_ieegBasics
addpath(genpath('path/to/eeglab/functions')); % https://github.com/sccn/eeglab/tree/develop/functions
addpath(genpath('path/to/colormaps'));

%% ~~~~~~~~~~~ Specify paths for saving the figures ~~~~~~~~~~~~~~~
SAVE_FIG = 0;
fig_path = '/path/to/figure_output_directory';

%% ~~~~~~~~~~~~~~~~~~~~~ Select DATA ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% Select files and specify filepaths
baseline = '/path/to/RANT_baselineS002R02.dat';
% Load impeance values if they exist
if isfile(replace(baseline,'.dat', '_info0.txt'))
    impedanceValues = extractImpedanceValues(replace(baseline,'.dat', '_info0.txt'));
else
    impedanceValues = [];
end

% -- Select BCI2000 BSEPs Data with GUI  ----
% !!! Select only recordings with the same parameters !!!!
data_folder = ['/Users/lampert.frederik/Documents/Mayo/CorTec/data' filesep];
[fnames, fpth] = uigetfile([data_folder '*.dat'],'Select BCI2000 .DAT FILE(s)','MultiSelect','on');
if isa(fnames,'char')   % prevents from getting error while loading data from single file
    fnames = {fnames};
end

%% ~~~~~~~~~~~~ SET PARAMETERS and Constants ~~~~~~~~~~~~~~~~~~~~~
% Select subject
sub = 'c05';

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

% Specify initial data curation parameters
bCh_thr = 3; % Bad channel detection threshold (in #std)
interpMethod = 'linear'; % Specify interpolation method to replace lost samples

% Select spatial re-referencing - SELECT ONE
% reref = 'Native';
% reref = 'CAR';         % Common average
reref = 'BP';
% reref = 'leadSpecific'; RR_chnls = [1 9 17 25]
% reref = 'leadCAR';   % Custom Common Average for cortec Grids

%  Filtering Parameters
notch_filt = false;    % Notch filtering
stopWidth = 5;
HP_filt = false;       % Band-pass filtering
HP_cut = 0.5;

% Set time limits for analysis
t_lims = [-0.5 1.5];

% Select option to detect stimulation onset:
%   1) 'ImplantStimulation': state vector, when stim was delivered via volatile commands or persistent function
%   2) 'Internal': When stimulation specified by persistent command, pause and repetitions
stimTrig = 'ImplantStimulation'; % 'Internal'


%% ~~~~~~~~~~~~~~~~~~~~~~~~~ Aggregate data ~~~~~~~~~~~~~~~~~~~~~~~~~~~~
signal = [];
StimulusCode = [];
ImplantStimulation = [];
PhaseInSequence = [];
ImplantLostSample = [];

stimPulses_ref = [];
stimTriggers_ref = [];
T_inter_stim_ref = [];
amp_f_ref = [];
FS_ref = [];
ref_ch_ref = [];
amp_f_options = {'57.5dB','51.5dB','45.5dB','39.5dB'};

for f = 1:numel(fnames)
    % ---------------- Load and extract data -------------------
    [fsignal, states, parameters, ftotal_samples] = ...
        LoadBCI2kData(fullfile(fpth, fnames{f}), interpMethod);

    signal = [signal; double(fsignal)];
    StimulusCode = [StimulusCode; double(states.StimulusCode)];
    ImplantStimulation = [ImplantStimulation; double(states.ImplantStimulation)];
    PhaseInSequence = [PhaseInSequence; double(states.PhaseInSequence)];
    ImplantLostSample = [ImplantLostSample; double(states.ImplantLostSample)];

    % Extract stimulation parameters
    if parameters.EnableStimulation.NumericValue == 1
        stimPulses = parseBCI2000_parameters_values(parameters, 'StimulationPulses');
        stimTriggers = parseBCI2000_parameters_values(parameters, 'StimulationTriggers');
    else
        stimPulses = table();
        stimTriggers = table();
    end

    % Determine interstimulus interval
    if parameters.ISIMinDuration.NumericValue ~= 0 && parameters.PostSequenceDuration.NumericValue == 0
        T_inter_stim = parameters.ISIMinDuration.NumericValue;
    elseif parameters.ISIMinDuration.NumericValue == 0 && parameters.PostSequenceDuration.NumericValue ~= 0
        T_inter_stim = parameters.PostSequenceDuration.NumericValue;
    else
        warning('Not able to determine the interstimulus interval from ISIMinDuration or PostSequenceDuration')
        T_inter_stim = NaN;
    end

    % Extract recording parameters
    amp_f = amp_f_options{parameters.AmplificationFactor.NumericValue + 1};
    FS = parameters.SamplingRate.NumericValue;
    ref_ch = parameters.ReferenceCh.NumericValue;

    % Compare parameters across files
    if f == 1
        stimPulses_ref = stimPulses;
        stimTriggers_ref = stimTriggers;
        T_inter_stim_ref = T_inter_stim;
        amp_f_ref = amp_f;
        FS_ref = FS;
        ref_ch_ref = ref_ch;
    else
        if ~isequaln(stimPulses, stimPulses_ref)
            error('StimulationPulses do not match across files. Mismatch found in file %s.', fnames{f});
        end
        if ~isequaln(stimTriggers, stimTriggers_ref)
            error('StimulationTriggers do not match across files. Mismatch found in file %s.', fnames{f});
        end
        if ~strcmp(amp_f, amp_f_ref)
            error('Amplification factor does not match across files. Mismatch found in file %s.', fnames{f});
        end
        if ~isequaln(T_inter_stim, T_inter_stim_ref)
            error('Interstimulus interval does not match across files. Mismatch found in file %s.', fnames{f});
        end
        if ref_ch ~= ref_ch_ref
            error('Reference channel does not match across files. Mismatch found in file %s.', fnames{f});
        end
    end
end

% Assign common parameters after aggregation
stimPulses = stimPulses_ref;
stimTriggers = stimTriggers_ref;
T_inter_stim = T_inter_stim_ref;
amp_f = amp_f_ref;
FS = FS_ref;
ref_ch = ref_ch_ref;
total_samples = size(signal,1);

clear stimPulses_ref stimTriggers_ref T_inter_stim_ref amp_f_ref FS_ref ref_ch_ref

%% ~~~~~~~~~~~~~~~~~ Pre-proccess the data ~~~~~~~~~~~~~~~~~~
BSEPs = struct(); % Initialize structre for data

% ---------- Bad chanel detection -------------
% Identify bad channels from pre-run sequence period
[bad_chnls, badChannelInfo] = identifyBadChannels(signal(PhaseInSequence==0,:), cnstrct, FS,3, knownBadChnls, impedanceValues, 1);

% ---------- Spatial re-referncing ----------
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

% ---------- Temporal filtering ----------
if  HP_filt % FIR high-pass filter using EEGLAB
    dataFilt = pop_eegfiltnew(proc_signal.', ... % EEGLAB expects data as [channels x samples]
        'hicutoff', HP_cut);
    proc_signal = dataFilt.';
end

%% ~~~~~~~~~~~~~~~~ Extract Stimulation periods (chop  data) ~~~~~~~~~~~~~~~~
% -------------------- Time windowing --------------------
tIdx = round(t_lims(1)*FS) : round(t_lims(2)*FS);

% Check what pulses were appied, which triggers were used
[segmentDur, triggerCode] = runlength(StimulusCode, total_samples);
tigerOnsets = [1; cumsum(segmentDur(1:end-1)) + 1];
stimOnsets = find(diff(ImplantStimulation) == 1);       % states.ImplantStimulation in uint8 (allowing only positive int)
stimCodes = double(StimulusCode(stimOnsets)); % Find stim onsets only for current StimulusCode and function IDS

switch stimTrig
    case 'ImplantStimulation'
        stimResponseWndw = [stimOnsets+tIdx(1), stimOnsets+tIdx(end)];
        validResp = stimResponseWndw(:,1) >= 1 & stimResponseWndw(:,2) <= total_samples;
        stimOnsets = stimOnsets(validResp);
        stimCodes = stimCodes(validResp);
        stimResponseWndw = stimResponseWndw(validResp,:);

        % Extract the stimulation response
        response = arrayfun(@(i) proc_signal(stimResponseWndw(i,1):stimResponseWndw(i,2),:), ...
            1:size(stimResponseWndw,1), 'UniformOutput', false); % get the respose within the window

    case 'Internal'
        response = {};
        resp_idx = [];
        respStimCode = [];

        for sIdx = 1:length(stimOnsets)
            trigger = stimCodes(sIdx);
            triggerExpr = stimTriggers{"Trigger", :};
            triggerValues = cellfun(@(x) sscanf(char(x), 'StimulusCode==%d'), triggerExpr);
            colIdx = find(triggerValues == trigger, 1);
            if isempty(colIdx)
                warning('No matching stimulation trigger found for StimulusCode == %d.', trigger);
                continue
            end

            funIDs = cell2mat(stimTriggers{"FunctionID_s_", colIdx});
            N = cell2mat(stimTriggers{"CommandRepetition", colIdx});

            stimT_us = 0;
            for fID = 1:numel(funIDs)
                pulseFuncIDs = cell2mat(stimPulses{"FunctionID", :});
                fIdx = find(pulseFuncIDs == funIDs(fID), 1);
                pulseDur = cell2mat(stimPulses{"Pulse duration_us_", fIdx});
                dz0 = cell2mat(stimPulses{"Dead zone 0_us_", fIdx});
                dz1 = cell2mat(stimPulses{"Dead zone 1_us_", fIdx});
                stimT_us = stimT_us + (5 * pulseDur) + (2 * dz0) + dz1; % Approximate full pulse-function duration in microseconds
            end
            stimT_samples = max(1, round((stimT_us / 1e6) * FS));

            for rep = 1:N
                thisStimOnset = stimOnsets(sIdx) + (rep-1) * stimT_samples;
                thisWndw = [thisStimOnset + tIdx(1), thisStimOnset + tIdx(end)];
                if thisWndw(1) < 1 || thisWndw(2) > total_samples
                    continue
                end
                response{end+1} = proc_signal(thisWndw(1):thisWndw(2), :); %#ok<SAGROW>
                resp_idx(end+1,:) = thisWndw; %#ok<SAGROW>
                respStimCode(end+1,1) = trigger; %#ok<SAGROW>
            end
        end

        stimCodes = respStimCode;
        stimResponseWndw = resp_idx;
end

%Extract unique stimulus codes applied
uniqueStimCodes = unique(stimCodes(:).');
uniqueStimCodes(isnan(uniqueStimCodes) | uniqueStimCodes == 0) = [];

for sc = 1:numel(uniqueStimCodes)
    triggerExpr = stimTriggers{"Trigger", :};
    scField = sprintf('StimulusCode%d',uniqueStimCodes(sc));
    triggerValues = cellfun(@(x) sscanf(char(x), 'StimulusCode==%d'), triggerExpr);
    trigCol = find(triggerValues == sc, 1);
    funIDs = cell2mat(stimTriggers{2, trigCol});
    stimChannels = struct();
    for ff = 1:numel(funIDs)
        pulseFuncIDs = stimPulses{"FunctionID", :};
        if iscell(pulseFuncIDs)
            pulseFuncIDs = cell2mat(pulseFuncIDs);
        end
        pulseCol = find(pulseFuncIDs == funIDs(ff), 1);
        anodes = stimPulses{6, pulseCol};
        cathodes = stimPulses{7, pulseCol};
        if iscell(anodes)
            stimChannels(ff).(scField).Anodes = cell2mat(anodes);
        else
            stimChannels(ff).(scField).Anodes = anodes;
        end
        if iscell(cathodes)
            stimChannels(ff).(scField).Cathodes = cell2mat(cathodes);
        else
            stimChannels(ff).(scField).Cathodes = cathodes;
        end
        
    end
end

%% ~~~~~~~~~~~~~~~~~~~~~~~~~~~ Plot response ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% Plotting parameters
alp = 0.2;                            % set alpha for plotting                                                                  % set offset (in ms) to display signal preceeding the stimulation

% 1st figure, raw signal and stm onsets
fig = figure('Position', [100 100 1500 800], 'Name','Stim Onsets');
plot(proc_signal, 'Color', [0 0 0 alp],'LineWidth',1);
hold on
stem(stimOnsets, repelem(max(signal,[],'all'),numel(stimOnsets)), 'Color','r')
hold off


c = ceil(sqrt(numel(ch_names)));  % calculate number of columns
r = ceil(numel(ch_names) / c);          % calculate number of rows
for sc = 1:numel(uniqueStimCodes)
    resp_data_cell = response(stimCodes==uniqueStimCodes(sc));
    resp_data_mat = reshape(cell2mat(resp_data_cell),[size(tIdx,2), numel(ch_names), numel(resp_data_cell)]); % convert resp_data from cell to matrix
    t_win = tIdx /FS;
    scField = sprintf('StimulusCode%d',uniqueStimCodes(sc));
    source_ch = stimChannels.(scField).Anodes;
    dest_ch =  stimChannels.(scField).Cathodes;

    % 2nd figure, average responses
    fig = figure('Position', [100 100 1500 800], 'Name', 'Average Response');
    tLayout = tiledlayout(r, c, 'TileSpacing', 'compact');
    for ch=1:size(resp_data_mat,2)
        nexttile
        plot(t_win, squeeze(resp_data_mat(:,ch,:)), 'Color',[0 0 0 alp]); hold on
        plot(t_win, mean(resp_data_mat(:,ch,:),3),'Color','red', 'LineWidth',1.5)
        xline(0,'--m'); % 'Stim onset'
        title(ch_names{ch});
        xlabel('Time [ms]'); ylabel('U [\muV]');
        axis tight;
    end

    srcStr = sprintf('%d ', source_ch);
    dstStr = sprintf('%d ', dest_ch);
    sgtitle(sprintf('Average response for %s: (Stim channels: [%s] → [%s])', ...
        scField, strtrim(srcStr), strtrim(dstStr)), 'FontSize',16)

    if SAVE_FIG == 1
        saveas(gcf,fullfile(fig_path, [fl_id '_' num2str(source_ch) '-' num2str(dest_ch) '.png']))
        % print(gcf,fullfile(fig_path, [fl_id '_' num2str(source_ch) '-' num2str(dest_ch)]),'-dpng','-r300')
        % print(gcf,fullfile(fig_path, [fl_id '_' num2str(source_ch) '-' num2str(dest_ch)]),'-dsvg','-r300','-vector')
    end

end
% save(replace(fl_pth,'.dat','_BSEPs.mat'), 'BSEPs');

clear plt sc srcStr dstStr

%% Pack everything into structure
% General Info
BSEPs.RecordingInfo.FileNames = fnames;
BSEPs.RecordingInfo.FilePath = fpth;
BSEPs.RecordingInfo.FS = FS;
BSEPs.RecordingInfo.AmplificationFactor = amp_f;
BSEPs.RecordingInfo.ReferenceChannel = ref_ch;
BSEPs.RecordingInfo.TotalSamples = total_samples;
BSEPs.RecordingInfo.ChannelNames = ch_names;
BSEPs.RecordingInfo.ImpedanceValues = impedanceValues;

% General Preprocessing
BSEPs.PreprocessingInfo.Reref = reref;
BSEPs.PreprocessingInfo.Construct = cnstrct;
BSEPs.PreprocessingInfo.BadChannels = bad_chnls;
BSEPs.PreprocessingInfo.BadChannelInfo = badChannelInfo;
BSEPs.PreprocessingInfo.HP_filt = HP_filt;
BSEPs.PreprocessingInfo.HP_cut = HP_cut;
BSEPs.PreprocessingInfo.Notch_filt = notch_filt;
% Data extraction parameters
BSEPs.PreprocessingInfo.ResponseWindowIndices = tIdx;
BSEPs.PreprocessingInfo.ResponseWindowLims = t_lims;

% Signals
BSEPs.Signal.Signal = proc_signal;
BSEPs.Signal.StimulusCode = StimulusCode;
BSEPs.Signal.ImplantStimulation = ImplantStimulation;
BSEPs.Signal.PhaseInSequence = PhaseInSequence;
BSEPs.Signal.ImplantLostSample =ImplantLostSample;

% Stim parameters
BSEPs.StimParameters.StimulationPulses = stimPulses;
BSEPs.StimParameters.StimulationTriggers = stimTriggers;
BSEPs.StimParameters.InterStimulusInterval = T_inter_stim;
BSEPs.StimParameters.StimTriggerMode = stimTrig;
BSEPs.StimParameters.StimChannels = stimChannels;

for sc = 1:numel(uniqueStimCodes)
    scField = sprintf('StimulusCode%d',uniqueStimCodes(sc));
    scMask = stimCodes == uniqueStimCodes(sc);
    BSEPs.Responses.(scField).Signal = response(scMask);
    BSEPs.Responses.(scField).Info.StimCode = sc;
    BSEPs.Responses.(scField).Info.StimOnsets = stimOnsets(scMask);
    BSEPs.Responses.(scField).Info.ResponseIdx = stimResponseWndw(scMask,:);

end

clearvars -except BSEPs FS fnames fpth fig_path SAVE_FIG reref baseline impedanceValues data_folder dog cnstrct spec knownBadChnls  t_lims  cm
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

parts = regexp(fnames{1}, 'S\d+', 'split');
fname = matlab.lang.makeValidName(parts{1}); 

if save_flag == true
    save(fullfile(fpth, [fname '_chopped_BSEPs']), 'BSEPs');
elseif save_flag == false
    disp('File will NOT be saved')
else
    error('BCI2000 dat2mat CNEL save: Could not save file %s',fname);
end


