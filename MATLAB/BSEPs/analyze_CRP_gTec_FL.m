% Application of CRP method on single-pulse stimuation data measured with
% gTec in canines
%%
clear all; close all; clc;
%% ~~~~~~~~~~~~ ADD PATHS ~~~~~~~~~~~~~~~
addpath('path/to/BCI2000/mex');
addpath('path/to/crp_scripts'); % https://github.com/kaijmiller/crp_scripts
addpath('path/to/mnl_ieegBasics/functions'); % https://github.com/MultimodalNeuroimagingLab/mnl_ieegBasics
addpath(genpath('path/to/eeglab/functions')); % https://github.com/sccn/eeglab/tree/develop/functions
addpath(genpath('path/to/colormaps'));
%% ~~~~~~~~~~~ SPECIFY PATHS ~~~~~~~~~~~~~
default_pth = '/path/to/data';
figpath = '/path/to/figure_output_directory';
%% ~~~~~~~~~~~ Parameters ~~~~~~~~~~~~~
figsave=0;
keepChnls = 1:24;
t_lims = [-0.5 2];
reref = 'DP'; %'DP'; % 'native'
cm = magma();
%% ~~~~~~~~~~~ Select File ~~~~~~~~~~~~~~~~~
[file, location] = uigetfile( ...
    {'*.dat;*.mat', 'Data files (*.dat, *.mat)'}, ...
    'Select a File');

if isequal(file,0)
    return; % user cancelled
end
[~, ~, ext] = fileparts(file);
%% ~~~~~~~~~~~ LOAD DATA ~~~~~~~~~~~~~~~~~
switch lower(ext)
    case '.dat' % If new data convert
        BSEPs = extractBSEPs_gtec(fullfile(location, file), [-0.5 2], keepChnls, reref);

        % Create output file name
        save_file = fullfile(location, ['chopped_' reref '_' strrep(file, '.dat', '.mat')]);

        if exist(save_file, 'file')
            % Ask user whether to overwrite
            choice = questdlg( ...
                sprintf('File already exists:\n\n%s\n\nOverwrite?', save_file), ...
                'Overwrite file?', ...
                'Yes', 'No', 'No');

            if ~strcmp(choice, 'Yes')
                disp('File not overwritten.');
                return;
            end
        end
        % Save converted data
        save(save_file, 'BSEPs', 't_lims', 'keepChnls', '-v7.3');
        
    case '.mat'
        load(fullfile(location, file)); % If preprocessed data load
    otherwise
        error('Unsupported file type: %s', ext);
end

clear choice
%% ~~~~~~~~ Visualize chopped siganal~~~~~~~~~~~
t = BSEPs.StimResponse.Time;
for i = 1:numel(BSEPs.StimResponse.Responses)
    respMat = BSEPs.StimResponse.Responses{i};
    
    figure('Name', sprintf('Stim Pair %d-%d', BSEPs.StimResponse.StimPair{i}), ...
        'Position', [100 100 1500 800]);
    tiledlayout('flow')
    for ch = 1:size(respMat,3)
        nexttile
        plot(t,squeeze(respMat(:,:,ch)),'Color',0.5*[1 1 1]);
        hold on
        plot(t,mean(respMat(:,:,ch),2),'k','LineWidth',1)
        title(BSEPs.Preprocessing.ChannelLabels(ch))
        xlim('tight'); ylim([-1e3 1e3])
        xlabel('Time (s)')
        ylabel('Amplitude (\muV)')
    end
end
clear t respMat i ch
%% ~~~~~~~~ MAIN LOOP : GET CRP PARAMETERS ~~~~~~~~~~~
if ~exist('BSEPs', 'var') || isempty(fieldnames(BSEPs))
    %Load Chopped data
    load(fullfile(location, file)); 
end

% ~~~~~~~~~ Setup parameters for CRP ~~~~~~~~~~~~~
t_win = BSEPs.StimResponse.Time;
fname = BSEPs.FileID; %21
opts.srate            = BSEPs.RecordingParameters.FS;
opts.stim_time        = 0;        % seconds
opts.t1               = 0.03;    % post-stim start (s)
opts.t2               = 1.0;      % post-stim end (s)
opts.t_baseline1      = 0.4;     % Start of the baseline (T_stim - t_baseline1*FS);
opts.t_baseline2      = 0.005;    % End of the baseline (T_stim - t_baseline2*FS);
opts.p_thresh         = 1e-4;     % conditions to be potentially artifactual
opts.artifact_removal = true;
% COMMENT ONE OPTION IN AND ONE OUT: if doing artifact removal, test for significance at response duration or at full datalength sent?
opts.art_interval     = 'tR';     % 'tR' or 'full' 

% ~~~~~~~~ Run CRP on the stim epochs ~~~~~~~~~~~
CRP_struct = struct;
CRP_struct.FileID = fname;
chNames = BSEPs.Preprocessing.ChannelLabels;
CRP_struct.ChannelNames = chNames;
CRP_struct.inputParameters = opts;
tic;
for i = 1:numel(BSEPs.StimResponse.Responses)
    V = BSEPs.StimResponse.Responses{i,1};

    if size(V,2) < 9 % get rid of stim pairs which were aborted before 9 stims
        return
    end

    stim_pair = BSEPs.StimResponse.StimPair{i};
    stimPair_str = sprintf('StimPair_%02d_%02d',stim_pair(1), stim_pair(2));
    fprintf('Processing %s\n\r',stimPair_str);
    for chan = 1:size(V,3)
        ch_nums = sscanf(chNames{chan}, 'Ch%d-Ch%d');
        if any(ch_nums == stim_pair(1) | ch_nums == stim_pair(2)) % Check if the analyzed channel does not contain Stim channel
            continue
        end
        [crp_projs,crp_parms,bad_trials,baselines]=CRP_call_dbs_FL(V(:,:,chan),t_win,opts);
        if strcmp(BSEPs.Preprocessing.ReRefType,'DP')
           hlp_key = strrep(chNames{chan},'-','_');
        else
           hlp_key = chNames{chan};
        end
        CRP_struct.(stimPair_str).(hlp_key).crp_parms = crp_parms;
        CRP_struct.(stimPair_str).(hlp_key).crp_projs = crp_projs;
        CRP_struct.(stimPair_str).(hlp_key).bad_trials = bad_trials;
        CRP_struct.(stimPair_str).(hlp_key).trial_baselines = baselines;
    end
end
clear hlp_key crp_parms crp_projs V baselines
toc;
save_file = strrep(save_file,'chopped','CRP');
save(save_file, 'CRP_struct', 't_lims', 'keepChnls', '-v7.3');

%% ~~~~~~~~~~~~~~~ Get the best responses ~~~~~~~~~~
if ~exist('CRP_struct', 'var') || isempty(fieldnames(CRP_struct))
    %Load CRP data
    load(strrep(fullfile(location,file),'chopped','CRP'));
end

% Extract parameters and sort them in ascending order
parameter = 'al_p';%'al_p';%'Vsnr'
stimLeadIdxs = 1:9;
results = {};  % Cell array to hold {StimPairName, ChannelName, PValue}
stimPairFields = fieldnames(CRP_struct);

for i = 1:numel(stimPairFields)
    stimPairName = stimPairFields{i};
    if ~startsWith(stimPairName,'StimPair') % Ommit any field that is not the stim field
        continue
    end
    stimPairStruct = CRP_struct.(stimPairName);
    
    chFields = fieldnames(stimPairStruct);
    
    for j = 1:numel(chFields)
        chName = chFields{j};
        chStruct = stimPairStruct.(chName);
        
        % Check if nested field exists
        if isfield(chStruct, 'crp_projs') && isfield(chStruct.crp_parms, parameter)
            parVal = chStruct.crp_parms.(parameter);
            % Average / median
            parVal = median(parVal);
            % Store: {StimPair, Channel, PValue}
            results(end+1, :) = {stimPairName, chName, parVal};
        end
    end
end

% Remove all the responses recorded on the stim lead
channelNums = cellfun(@(s) sscanf(s, 'Ch%d_Ch%d').', results(:,2), 'UniformOutput', false);
toRemove = cellfun(@(v) any(ismember(v,stimLeadIdxs)), channelNums);
results = results(~toRemove, :);

% Convert to table for easy sorting
T = cell2table(results, 'VariableNames', {'StimPair', 'Channel', parameter});

% Sort by parameter value
T = sortrows(T, parameter, 'descend');

% clear stimPairStruct chStruct stimPairName parVal channelNums toRemove

% ~~~~~~~~~~~~~~~~~ Plot heatmap for on alpha prime, SNR, etc ~~~~~~~~~~~~
figure
% stim pair vs channel
h = heatmap(T,'Channel','StimPair','ColorVariable',parameter, 'ColorMethod','median');
% Replace underscores in display labels (X and Y axes)
h.XDisplayLabels = strrep(h.XDisplayLabels, '_', '-');
h.YDisplayLabels = strrep(h.YDisplayLabels, '_', '-');
h.Colormap = cm;

%% ~~~~~~~~~~~~~~~~~ Plot the CRP parametrization ~~~~~~~~~~~~~~
% ================= User parameters =================
N = 50;                  % number of responses to plot
thresh_pct = 0.98;      % this estimate of uncertainty in \tau_R is chosen for when projection magnitude exceeds 98% of maximum.

% ===================================================

if ~exist('t_win')
    t_win = BSEPs.StimResponse.Time;
end

% CRP time window
t1 = CRP_struct.inputParameters.t1;
t2 = CRP_struct.inputParameters.t2;
tpts = t_win(t_win > t1 & t_win <= t2);

for iResp = 1:min(N, height(T))
    % -------- Extract identifiers from table --------
    stimPair_str = T.StimPair{iResp};
    % stimPair_str = 'StimPair_07_04';
    ch_str       = T.Channel{iResp};
    % ch_str       = 'Ch19_Ch20';

    stimPair = sscanf(stimPair_str, 'StimPair_%d_%d')';

    % -------- Resolve indices --------
    stimPair_idx = find(cellfun(@(x) isequal(x, stimPair), ...
        BSEPs.StimResponse.StimPair));

    ch_label = strrep(ch_str, '_', '-');
    ch_idx   = find(strcmp(BSEPs.Preprocessing.ChannelLabels, ch_label));

    % -------- Load raw data --------
    V_raw = BSEPs.StimResponse.Responses{stimPair_idx}(:,:,ch_idx);
    crp   = CRP_struct.(stimPair_str).(ch_str);

    baselines = crp.trial_baselines;

    % Remove bad trials
    if ~isempty(crp.bad_trials)
        V_raw(:, crp.bad_trials) = [];
        baselines(crp.bad_trials) = [];
    end

    % Baseline subtraction
    V_raw = V_raw - baselines(:)';

    % -------- Uncertainty bounds on τR --------
    avg_proj = crp.crp_projs.mean_proj_profile;
    [mp, c]  = max(avg_proj);

    % lower bound
    b = 1 + find( ...
        avg_proj(2:end) > thresh_pct * mp & ...
        avg_proj(1:end-1) <= thresh_pct * mp );
    b(b > c) = [];
    if isempty(b)
        lb = 1;
    else
        lb = b(end);
    end
    

    % upper bound
    b = find( ...
        avg_proj(1:end-1) > thresh_pct * mp & ...
        avg_proj(2:end) <= thresh_pct * mp );
    b(b < c) = [];

    if isempty(b)
        hb = numel(avg_proj);
    else
        hb = b(1);
    end

    err_low_t  = crp.crp_projs.proj_tpts(lb);
    err_high_t = crp.crp_projs.proj_tpts(hb);

    % -------- Plotting --------
    figure('Name', sprintf('%s | %s', stimPair_str, ch_str), ...
           'NumberTitle', 'off');

    % Overall main title
    sgtitle(sprintf('Data traces and time-resolved projection magnitudes \n %s',... 
    [strrep(stimPair_str,'_','-'), '; ' strrep(ch_str,'_','-')]));

    % ---- Raw traces ----
    subplot(1,2,1); hold on
    plot(t_win, V_raw, 'Color', 0.5*[1 1 1])
    % ---- Plot baseline ----
    plot([t_win(1) t_win(end)], [0 0], 'Color', 0.3*[1 1 1])
    % ---- Plot/Highlight Response Time tR on average trace ----
    plot(crp.crp_parms.parms_times, ...
         crp.crp_parms.avg_trace_tR, ...
         'Color', [1 1 0], 'LineWidth', 4)
    % ---- Plot average trace ----
    plot(tpts, crp.crp_projs.avg_trace_input, 'k', 'LineWidth', 2)

    title(sprintf('t = %.2f, p = %.2e', ...
        crp.crp_projs.t_value_tR, crp.crp_projs.p_value_tR))
    % ---- Clip only to the significant response ---
    xlim([tpts(1) tpts(end)])
    % xlim([tpts(1) tpts(find(tpts <= 0.25, 1, 'last'))])
    % xlim([tpts(1) tpts(find(tpts <= crp.crp_parms.tR, 1, 'last'))])
    xlabel('Time (s)')
    ylabel('Voltage (\muV)')
    box off

    % ---- Projection magnitudes ----
    subplot(1,2,2); hold on
    plot(crp.crp_projs.proj_tpts, crp.crp_projs.S_all, ...
        'Color', 0.3*[1 1 1], 'LineWidth', 0.5)

    plot(crp.crp_projs.proj_tpts, avg_proj, 'k', 'LineWidth', 2)
    plot(crp.crp_parms.tR, max(avg_proj), 'ro', 'LineWidth', 2)

    a = diff(ylim) * 0.03;
    plot(err_low_t  * [1 1], max(avg_proj) + a*[-1 1], 'r', 'LineWidth', 2)
    plot(err_high_t * [1 1], max(avg_proj) + a*[-1 1], 'r', 'LineWidth', 2)
    xlim([tpts(1) tpts(find(tpts <= 0.25, 1, 'last'))])
    % xlim([tpts(1) tpts(find(tpts <= crp.crp_parms.tR, 1, 'last'))])
    xlabel('Time (s)')
    ylabel('Proj Mag (\muV s^{1/2})')
    title(sprintf('mean S = %.2f', max(avg_proj)))
    box off

end
%% 8 - Plot projection magnitudes of full epoch input and at response duration (\tau_R)
figure('Name','Cross-projection magnitudes at response duration and for full input','NumberTitle','off')
    
    subplot(1,2,1)
    plot(CRP_struct.(stimPair_str).(ch).crp_projs.S_all(:,end),'.','MarkerSize',15)
    ylabel('Proj Mag (\muV s^{1/2})'),xlabel('sample number'), box off
    title('Projection magnitude of full interval input')
    
    subplot(1,2,2)   
    plot(CRP_struct.(stimPair_str).(ch).crp_projs.S_all(:,find(CRP_struct.(stimPair_str).(ch).crp_projs.proj_tpts==CRP_struct.(stimPair_str).(ch).crp_parms.tR)),'.','MarkerSize',15)
    ylabel('Proj Mag (\muV s^{1/2})'),xlabel('sample number'), box off
    title('Projection magnitude at response duration \tau_R')
    
    if figsave=='y'
        kjm_printfig([figpath filesep fname '_projmags_tR_and_input'],[10  6])
    end    
    
    

%% display parameters
    disp('- - parameterizations - -')
    disp(['mean projection (S) at response duration: ' num2str(max(crp_projs.mean_proj_profile))])
    disp(['response duration, \tau_R: ' num2str(crp_parms.tR)])
    disp(['mean alpha: ' num2str(mean(crp_parms.al))])
    disp(['mean alpha-prime: ' num2str(mean(crp_parms.al)/(length(crp_parms.C).^.5))])
    disp(['mean epep_root: ' num2str(mean(crp_parms.epep_root))])
    disp(['mean SNR: ' num2str(mean(crp_parms.Vsnr))])
    disp(['mean explained variance: ' num2str(mean(crp_parms.expl_var))])
    disp(['Extraction significance t=' num2str(crp_projs.t_value_tR) ', p=' num2str(crp_projs.p_value_tR)])
    [~,tmp_p]=ttest(crp_parms.al);
    disp(['alpha coefficient significance p=' num2str(tmp_p)]), clear tmp_p
    disp('- - - - - - - - - - - - - - - - - - - - - - - - - - - -')

