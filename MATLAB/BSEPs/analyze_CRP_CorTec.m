% Script to run application of CRP method on single-pulse stimuation data measured with
% CorTec in canines

% Dependencies:
%    - .mat file obtained by runnin "extractBSEPs_CorTec"
%    - crp_scripts
%    - colormaps

%%
clear all; close all; clc;
%% ~~~~~~~~~~~~ ADD PATHS ~~~~~~~~~~~~~~~
addpath('path/to/BCI2000/mex');
addpath('path/to/crp_scripts'); % https://github.com/kaijmiller/crp_scripts
addpath(genpath('path/to/colormaps'));

%% ~~~~~~~~~~~ Specify paths for saving the figures ~~~~~~~~~~~~~~~
SAVE_FIG = 0;
figpath  = '/path/to/figure_output_directory';

%% ~~~~~~~~~~~ LOAD DATA ~~~~~~~~~~~~~~~~~
% Select files and specify filepaths
data_folder = ['/path/to/data' filesep];
% -- Select BSEPs Data with GUI  ----
[fname, fpth] = uigetfile([data_folder '*.mat'],'Select processed .mat FILE','MultiSelect','off');
load([fpth filesep fname]);

%% ~~~~~~~~~ Extract data ~~~~~~~~~~
FS = BSEPs.RecordingInfo.FS;
resp_wndw = BSEPs.PreprocessingInfo.ResponseWindowIndices;
t_win =  resp_wndw/FS;
ch_names = BSEPs.RecordingInfo.ChannelNames;
stimCodes = fieldnames(BSEPs.Responses);
stimChannels = BSEPs.StimParameters.StimChannels;
parts = regexp(BSEPs.RecordingInfo.FileNames{1}, 'S\d+', 'split');
figFname = matlab.lang.makeValidName(parts{1}); 

%% ~~~~~~~~ Visualize chopped siganal~~~~~~~~~~~

% Plotting parameters
alp = 0.2;                            % set alpha for plotting
c = ceil(sqrt(numel(ch_names)));  % calculate number of columns
r = ceil(numel(ch_names) / c);          % calculate number of rows

for sc = 1:numel(stimCodes)
    stimID = stimCodes{sc};
    resp_data_cell = BSEPs.Responses.(stimID).Signal;
    resp_data_mat = reshape(cell2mat(resp_data_cell),[size(resp_wndw,2), numel(ch_names), numel(resp_data_cell)]); % convert resp_data from cell to matrix
    source_ch = stimChannels.(stimID).Anodes;
    dest_ch = stimChannels.(stimID).Cathodes;

    % figure, average responses
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
        stimID, strtrim(srcStr), strtrim(dstStr)), 'FontSize',16)

    if SAVE_FIG
        print(gcf,fullfile(figpath, [figFname '_averageResponse']),'-dpng','-r300')
        % print(gcf,fullfile(figpath, [figFname 'averageResponse']),'-dsvg','-r300','-vector')
    end
end

clear alp c ch dest_ch dstStr r resp_data_cell resp_data_mat source_ch srcStr 

%% ~~~~~~~~ MAIN LOOP : GET CRP PARAMETERS ~~~~~~~~~~~
if ~exist('BSEPs', 'var') || isempty(fieldnames(BSEPs))
    %Load Chopped data
    load([fpth filesep fname]); 
end

% ~~~~~~~~~ Setup parameters for CRP ~~~~~~~~~~~~~
FS = BSEPs.RecordingInfo.FS;
t_win = BSEPs.PreprocessingInfo.ResponseWindowIndices/FS;
parts = regexp(BSEPs.RecordingInfo.FileNames{1}, 'S\d+', 'split');
fname = matlab.lang.makeValidName(parts{1}); 
opts.srate            = FS;
opts.stim_time        = 0;        % seconds
opts.t1               = 0.03;    % post-stim start (s)
opts.t2               = 1.5;      % post-stim end (s)
opts.t_baseline1      = 0.25;     % Start of the baseline (T_stim - t_baseline1*FS);
opts.t_baseline2      = 0.005;    % End of the baseline (T_stim - t_baseline2*FS);
opts.p_thresh         = 1e-4;     % conditions to be potentially artifactual
opts.artifact_removal = true;
% COMMENT ONE OPTION IN AND ONE OUT: if doing artifact removal, test for significance at response duration or at full datalength sent?
opts.art_interval     = 'tR';     % 'tR' or 'full' 

% ~~~~~~~~ Run CRP on the stim epochs ~~~~~~~~~~~
stimCodes = fieldnames(BSEPs.Responses);
CRP_struct = struct;
CRP_struct.FileID = fname;
chNames = BSEPs.RecordingInfo.ChannelNames;
CRP_struct.ChannelNames = chNames;
CRP_struct.inputParameters = opts;
tic;
for sc = 1:numel(stimCodes)
    stimID = stimCodes{sc};
    resp_data_cell = BSEPs.Responses.(stimID).Signal;

    V = reshape(cell2mat(resp_data_cell),[size(resp_wndw,2), numel(chNames), numel(resp_data_cell)]); % convert resp_data from cell to matrix
    if size(V,3) < 10 % get rid of stim pairs which were aborted before 9 stims
        return
    end

    source_ch = stimChannels.(stimID).Anodes;
    dest_ch =  stimChannels.(stimID).Cathodes;
    % stim_pair = BSEPs.StimResponse.StimPair{i};
    srcStr = sprintf(' %d', source_ch);
    dstStr = sprintf(' %d', dest_ch);
    fprintf('Processing response for %s: (Stim channels: [%s] → [%s])\n\r',stimID, srcStr, dstStr);

    for chan = 1:size(V,2)
        ch_nums = sscanf(chNames{chan}, 'Ch%d-%d');
        if any(ismember(ch_nums, [source_ch(:); dest_ch(:)])) % Check if the analyzed channel does not contain Stim channel
            continue
        end
        crp_data = squeeze(V(:,chan,:));
        [crp_projs,crp_parms,bad_trials,baselines]=CRP_call_dbs_FL(crp_data,t_win,opts);
        if strcmp(BSEPs.PreprocessingInfo.Reref,'BP')
            hlp_key = strrep(chNames{chan},'-','_');
        else
            hlp_key = chNames{chan};
        end
        CRP_struct.(stimID).(hlp_key).crp_parms = crp_parms;
        CRP_struct.(stimID).(hlp_key).crp_projs = crp_projs;
        CRP_struct.(stimID).(hlp_key).bad_trials = bad_trials;
        CRP_struct.(stimID).(hlp_key).trial_baselines = baselines;
    end

end
clear hlp_key crp_parms crp_projs V baselines
toc;
save_file = fullfile(fpth, [fname '_CRP']);
save(save_file, 'CRP_struct', 'resp_wndw', '-v7.3');


%% ~~~~~~~~~~~~~~~ Get the best responses ~~~~~~~~~~
if ~exist('CRP_struct', 'var') || isempty(fieldnames(CRP_struct))
    %Load CRP data
    load(fullfile(fpth,strrep(fname,'chopped_BSEPs','CRP')));
end

% Extract parameters and sort them in ascending order
parameter = 'Vsnr';%'al_p';%'Vsnr'
stimLeadIdxs = [];
results = {};  % Cell array to hold {StimPairName, ChannelName, PValue}
stimPairFields = fieldnames(CRP_struct);
% Set colormap
cm = magma();

for i = 1:numel(stimPairFields)
    stimID = stimPairFields{i};
    if ~startsWith(stimID,'StimulusCode') % Ommit any field that is not the stim field
        continue
    end
    stimIdStruct = CRP_struct.(stimID);
    chFields = fieldnames(stimIdStruct);
    
    for j = 1:numel(chFields)
        chName = chFields{j};
        chStruct = stimIdStruct.(chName);
        
        % Check if nested field exists
        if isfield(chStruct, 'crp_projs') && isfield(chStruct.crp_parms, parameter)
            parVal = chStruct.crp_parms.(parameter);
            % Average / median
            parVal = median(parVal);
            % Store: {StimPair, Channel, PValue}
            results(end+1, :) = {stimID, chName, parVal};
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

% Housekeeping
clear  parameter stimLeadIdxs  toRemove stimPairFields  stimIdStruct chFields chStruct parVal
%% ~~~~~~~~~~~~~~~~~ Plot the CRP parametrization ~~~~~~~~~~~~~~
% ================= User parameters =================
N = 20;                  % number of responses to plot
thresh_pct = 0.98;      % this estimate of uncertainty in \tau_R is chosen for when projection magnitude exceeds 98% of maximum.

% ===================================================
if ~exist('t_win')
    t_win = BSEPs.PreprocessingInfo.ResponseWindowIndices;
end
chNames = BSEPs.RecordingInfo.ChannelNames;

% CRP time window
t1 = CRP_struct.inputParameters.t1;
t2 = CRP_struct.inputParameters.t2;
tpts = t_win(t_win > t1 & t_win <= t2);

for iResp = 1:min(N, height(T))
    % -------- Extract identifiers from table --------
    stimID = T.StimPair{iResp};
    stimID = 'StimulusCode1';
    ch_str       = T.Channel{iResp};
    ch_label = strrep(ch_str, '_', '-');
    ch_idx   = find(strcmp(chNames, ch_label));
    source_ch = stimChannels.(stimID).Anodes;
    dest_ch =  stimChannels.(stimID).Cathodes;
    srcStr = sprintf('%d ', source_ch);
    dstStr = sprintf('%d ', dest_ch);

    % -------- Load raw data --------
    resp_data_cell = BSEPs.Responses.(stimID).Signal;
    V_raw = cell2mat(cellfun(@(x) x(:, ch_idx), resp_data_cell, 'UniformOutput', false));
    crp   = CRP_struct.(stimID).(ch_str);
 
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
    figure('Position', [100 100 1500 800],'Name', sprintf('%s | %s', stimID, ch_str), ...
           'NumberTitle', 'off');

    % Overall main title
    sgtitle(sprintf('Data traces and time-resolved projection magnitudes \n %s',... 
    [strrep(stimID,'_','-'), '; ' strrep(ch_str,'_','-')]));

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

    % ylim([-3000 3000])

    title(sprintf('t = %.2f, p = %.2e', ...
        crp.crp_projs.t_value_tR, crp.crp_projs.p_value_tR))
    % ---- Clip only to the significant response ---
    % xlim([tpts(1) tpts(find(tpts <= 0.25, 1, 'last'))])
    % xlim([tpts(1) tpts(find(tpts <= crp.crp_parms.tR, 1, 'last'))])
    xlim([t_win(find(t_win==0)) t_win(end)])
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
    % xlim([tpts(1) tpts(find(tpts <= 0.25, 1, 'last'))])
    xlim([tpts(1) tpts(find(tpts <= crp.crp_parms.tR, 1, 'last'))])
    xlabel('Time (s)')
    ylabel('Proj Mag (\muV s^{1/2})')
    title(sprintf('mean S = %.2f', max(avg_proj)))
    box off

    % display parameters
    disp('- - parameterizations - -')
    fprintf('Response for %s: (Stim channels: [%s] → [%s])\n\r', strtrim(srcStr), strtrim(dstStr))
    disp(['mean projection (S) at response duration: ' num2str(max(crp.crp_projs.mean_proj_profile))])
    disp(['response duration, \tau_R: ' num2str(crp.crp_parms.tR)])
    disp(['mean alpha: ' num2str(mean(crp.crp_parms.al))])
    disp(['mean alpha-prime: ' num2str(mean(crp.crp_parms.al)/(length(crp.crp_parms.C).^.5))])
    disp(['mean epep_root: ' num2str(mean(crp.crp_parms.epep_root))])
    disp(['mean SNR: ' num2str(mean(crp.crp_parms.Vsnr))])
    disp(['mean explained variance: ' num2str(mean(crp.crp_parms.expl_var))])
    disp(['Extraction significance t=' num2str(crp.crp_projs.t_value_tR) ', p=' num2str(crp.crp_projs.p_value_tR)])
    [~,tmp_p]=ttest(crp.crp_parms.al);
    disp(['alpha coefficient significance p=' num2str(tmp_p)]), clear tmp_p
    disp('- - - - - - - - - - - - - - - - - - - - - - - - - - - -')

    if SAVE_FIG
        % print(gcf,fullfile(figpath, [CRP_struct.FileID   '_' stimID '_' ch_str]),'-dpng','-r300')
        print(gcf,fullfile(figpath, [CRP_struct.FileID   '_' stimID '_' ch_str]),'-dsvg','-r300','-vector')
    end
end


clear a b bad_trials baselines c ch ch_idx ch_label ch_names ch_nums ch_Str chan chNamedest_ch dstStr err_high_t err_low_t hb i iResp j k lb mp N sc source_ch

%% Plot Spectral map
win = hann(0.1 * FS);
nOverlap = 0.5*length(win);
freqSpec = 256;
smoothF = [5,5]; % Smoothing factor  across frequency and time

    figure('Position', [200 100 1000 800],'Name', sprintf('%s | %s', stimID, ch_str), ...
           'NumberTitle', 'off');

    % Overall main title
    sgtitle(sprintf('Data traces and TF map \n %s',... 
    [strrep(stimID,'_','-'), '; ' strrep(ch_str,'_','-')]));

    % ---- Raw traces ----
    subplot(2,1,1); hold on
    plot(t_win, V_raw, 'Color', 0.5*[1 1 1]) 

    % ---- Plot baseline ----
    plot([t_win(1) t_win(end)], [0 0], 'Color', 0.3*[1 1 1])
    % ---- Plot/Highlight Response Time tR on average trace ----
    plot(crp.crp_parms.parms_times, ...
         crp.crp_parms.avg_trace_tR, ...
         'Color', [1 1 0], 'LineWidth', 4)
    % ---- Plot average trace ----
    plot(tpts, crp.crp_projs.avg_trace_input, 'k', 'LineWidth', 2)

     % ---- TF Map ----
     subplot(2,1,2); hold on
     [S,F,T] = spectrogram(mean(V_raw,2), win, nOverlap, freqSpec, FS);
     
     % ------------------
     T = T + t_win(1);
     P = pow2db(abs(S).^2);
     % Smooth spectrogram
     P = smoothdata(P, 1, 'gaussian', smoothF(1));   % Smooth across frequency
     P = smoothdata(P, 2, 'gaussian', smoothF(2));   % Smooth across time
     imagesc(T, F, P);
     axis xy
     xline(0,'LineWidth',2,'Color','r')
     xlim([t_win(1) t_win(end)])
     xticks(t_win(1):0.5:t_win(end))
     xlabel('Time (s)')
     ylabel('Frequency (Hz)')
     % colormap(turbo)
     cb = colorbar('southoutside');
     cb.Label.String = 'Power (dB/Hz)';
     axis tight