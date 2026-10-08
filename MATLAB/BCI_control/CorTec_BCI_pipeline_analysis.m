clear all; close all; clc
%% Add Paths
addpath('./functions/');
addpath('.../violinplot'); % https://github.com/bastibe/Violinplot-Matlab.git
addpath('.../legendUnq'); % https://www.mathworks.com/matlabcentral/fileexchange/67646-legendunq
addpath('path/to/mnl_ieegBasics/functions'); % https://github.com/MultimodalNeuroimagingLab/mnl_ieegBasics

%% Load KJM colormaps
rwb = load('colormaps/rwb_colormap.mat');
rgb = load('colormaps/RYGyCB.mat');

%% Load the processed data
[fname, fpth] = uigetfile(['/path/to/data' filesep '*.mat'], ...
    'Select the INPUT .mat FILE','MultiSelect','off');      % UI selectfile to open
load([fpth filesep fname]);

%% Specify path where to save figures
SAVE_FIG = 0;
fig_pth = '/path/to/figure_output_directory';

%% Retrieve subject ID and Runs
subID = data.subjectID;
fld_names = fieldnames(data);
fld_names = fld_names(~ismember(fld_names, {'baseline', 'subjectID'}));
%% Plot whole BCI runs 
i = 2; % select index of the recording you want to plot
SR = data.(fld_names{i}).SamplingRate;
signal = data.(fld_names{i}).RerefSignal.Signal;
trialSpectFeatures = data.(fld_names{i}).RerefSignal.SpectralFeatures;
ch_names = data.(fld_names{i}).RerefSignal.Info.InputLabels;
BP_label = strjoin(ch_names, '-');
% StimCode 3 is down, StimCode 4 is UP
h = plotTrialSpectralBrowser(signal, trialSpectFeatures, SR, BP_label);
h.axTF.CLim = [-10 10];

%% Extract BCI accuracy and performance
normOutputIdx = 1; % Select which Normalzier channel was used to control BCI task 
targets = [];
performance = [];
nTrials = [];
avgProcessingOutput = [];
medProcessingOutput = [];

for i = 1:length(fld_names)
    BCI_task = data.(fld_names{i}).BCITask; % Extract the BCI task results

    % Get the results of the trials and their indices
    targets = [targets; BCI_task.targetCode];
    performance = [performance; BCI_task.performance];
    nTrials = [nTrials; BCI_task.nTrials];

    % Get the mean output of the porcessing for selected trial
    outField = sprintf('SGNL_NormalizerOutput_%d',normOutputIdx);
    normOutput = data.(fld_names{i}).BCI2000ProcessingOutputs.(outField);

    % Get the start and stop indices of the trials
    startIdx = BCI_task.startIdx;
    endIdx = BCI_task.endIdx;

    % Extract Normalizer values
    for tr = 1:numel(startIdx)
        avgProcessingOutput = [avgProcessingOutput; mean(normOutput(startIdx(tr):endIdx(tr)))];
        medProcessingOutput = [medProcessingOutput; median(normOutput(startIdx(tr):endIdx(tr)))];
    end

end
clear BCI_task tr outField normOutput startIdx endIdx


%% Plot BCI performance
% Select which statistic will be plotted
plotData = avgProcessingOutput; % medProcessingOutput

% Create figure
hFig = figure('Visible','on');
set(hFig, 'Position', [200, 100, 1400, 400]);
tl = tiledlayout(hFig, 1, 4, 'TileSpacing', 'compact', ...
    'Padding', 'compact');
lw = 1.5;
markerSize = 10;

% Main BCI trial-wise plot
ax1 = nexttile(tl, [1 3]);
hold(ax1, 'on');
yline(ax1, 0, 'k', 'HandleVisibility', 'off');
for tr = 1:length(plotData)
    if performance(tr) == 0 && targets(tr) == 3 % Unsuccessful lower
        plot(ax1, tr, plotData(tr), 'b*', 'LineWidth', lw, ...
            'MarkerSize', markerSize, 'DisplayName', 'Unsuccessful lower');
    elseif performance(tr) == 0 && targets(tr) == 4 % Unsuccessful upper
        plot(ax1, tr, plotData(tr), 'r*', 'LineWidth', lw, ...
            'MarkerSize', markerSize, 'DisplayName', 'Unsuccessful upper');
    elseif performance(tr) == 1 && targets(tr) == 3 % Successful lower
        plot(ax1, tr, plotData(tr), 'bo', 'LineWidth', lw, ...
            'MarkerSize', markerSize, 'DisplayName', 'Successful lower');
    elseif performance(tr) == 1 && targets(tr) == 4 % Successful upper
        plot(ax1, tr, plotData(tr), 'ro', 'LineWidth', lw, ...
            'MarkerSize', markerSize, 'DisplayName', 'Successful upper');
    elseif performance(tr) == -1 && targets(tr) == 3 % Incorrect lower
        plot(ax1, tr, plotData(tr), 'bx', 'LineWidth', lw, ...
            'MarkerSize', markerSize, 'DisplayName', 'Incorrect lower');
    elseif performance(tr) == -1 && targets(tr) == 4  % Incorrect upper
        plot(ax1, tr, plotData(tr), 'rx', 'LineWidth', lw, ...
            'MarkerSize', markerSize, 'DisplayName', 'Incorrect upper');
    else % fallback
        continue
    end
end

xline(ax1, cumsum(nTrials)+0.5, 'LineWidth', 0.75, ...
    'LineStyle', '--', 'Color', 'k', 'DisplayName', 'End of run');
xlim([0 sum(nTrials)+1])
ylabel(ax1, 'Mean Normalized Signal Processing Output');
xlabel(ax1, 'Trials');
title(ax1, 'BCI performance');
grid(ax1, 'on');
legend(ax1, legendUnq(hFig));

% Violin plot
ax2 = nexttile(tl, 4);
hold(ax2, 'on');

% Violin plot
hitsUp  = plotData((performance == 1) & (targets == 4));
hitsLow = plotData((performance == 1) & (targets == 3));
cats = [repelem({'Upper Hit'}, numel(hitsUp)), ...
        repelem({'Lower Hit'}, numel(hitsLow))];
violinplot([hitsUp(:); hitsLow(:)], cats);
ylabel(ax2, 'Normalized Output');
title(ax2, 'Hit distribution');
grid(ax2, 'on');

fontsize(hFig, 16, 'points');


%% Second violinplot with all the data
markerSize = 30;

hitsUp  = plotData((performance == 1) & (targets == 4));
hitsLow = plotData((performance == 1) & (targets == 3));
timeOutUp = plotData((performance == 0) & (targets == 4));
timeOutLow = plotData((performance == 0) & (targets == 3));
missUp = plotData((performance == -1) & (targets == 4));
missLow =  plotData((performance == -1) & (targets == 3));

allData = [hitsUp(:); hitsLow(:); ...
           missUp(:); missLow(:); ...
           timeOutUp(:); timeOutLow(:)];

cats = [repelem({'Correct Up'}, numel(hitsUp)), ...
        repelem({'Correct Low'}, numel(hitsLow)), ...
        repelem({'Miss Up'}, numel(missUp)), ...
        repelem({'Miss Low'}, numel(missLow)), ...
        repelem({'Timeout Up'}, numel(timeOutUp)), ...
        repelem({'Timeout Low'}, numel(timeOutLow))];

figure('Position', [200, 100, 1400, 400]);
hold 'on';
vs = violinplot([hitsUp(:); hitsLow(:); ...
    missUp(:); missLow(:); ...
    timeOutUp(:); timeOutLow(:)], ...
    cats, ...
    'MarkerSize', markerSize);
vs(3).ScatterPlot.Marker = 'x';
vs(4).ScatterPlot.Marker = '*';
ylabel('Normalized Output');
title('Hit distribution');
grid('on'); 
fontsize(16, 'points');

%% Extract and Plot overall trial PSDs
ch = 1; % Channel to plot
PSDs = [];
notchFreq = [60 120 180];
LCFreqs = [72.5 112.5];

yl = [0.01 10];
for i = 1:length(fld_names)
    trialSpectFeatures = data.(fld_names{i}).RerefSignal.SpectralFeatures;
    tsc = data.(fld_names{i}).RerefSignal.SpectralFeatures.trialStimCode;
    PSDs = cat(3,PSDs,trialSpectFeatures.trialPSD(:,:,tsc~=0));
end

freqBins = trialSpectFeatures.freqBins; % Extract Frequncies from last run
deltaF = min(diff(freqBins));
notchMask = find(any(freqBins<= notchFreq+deltaF &  freqBins>= notchFreq-deltaF,2)); % Find all the freq bins to remove
PSDs(notchMask,:,:) = nan;

hFig = figure('Position',[200, 100, 1400, 650]);
tl = tiledlayout(hFig, 1, 2, 'TileSpacing', 'compact', ...
    'Padding', 'compact');

% Left tile (PSDs across all trials)
ax1 = nexttile(tl);
semilogy(ax1, freqBins, mean(PSDs(:,ch,targets==3),3,'omitnan'), ...
    'Color','b', 'LineWidth',1.5, 'DisplayName','Avg PSD for low'), hold on
semilogy(ax1, freqBins, mean(PSDs(:,ch,targets==4),3,'omitnan'), ...
    'Color','r', 'LineWidth',1.5, 'DisplayName','Avg PSD for upper')
% Patch the broadband range
% yl = get(ax3, 'Ylim');
patch(ax1, [LCFreqs(1) LCFreqs(2) LCFreqs(2) LCFreqs(1)], ...
    [yl(1) yl(1) yl(2) yl(2)], ...
    [0.5 0.5 0.5], 'FaceAlpha', 0.18, 'EdgeColor', 'none', 'DisplayName', 'Broadband range');
xlabel(ax1, 'Frequency [Hz]'); 
ylabel(ax1, 'Power \muV^2 / Hz');
title(ax1, 'Average PSDs across all trials')
grid(ax1, 'on');
xlim([0 350])

% Right tile (PSDs of only successfull trials)
ax2 = nexttile(tl);
semilogy(ax2, freqBins, ...
    mean(PSDs(:,ch,(performance == 1) & (targets == 3)),3,'omitnan'), ...
    'Color','b', 'LineWidth',1.5, 'DisplayName','Avg PSD for low'), hold on
semilogy(ax2, freqBins, ...
    mean(PSDs(:,ch,(performance == 1) & (targets == 4)),3,'omitnan'), ...
    'Color','r', 'LineWidth',1.5, 'DisplayName','Avg PSD for upper')
% Patch the broadband range
patch(ax2, [LCFreqs(1) LCFreqs(2) LCFreqs(2) LCFreqs(1)], ...
    [yl(1) yl(1) yl(2) yl(2)], ...
    [0.5 0.5 0.5], 'FaceAlpha', 0.18, 'EdgeColor', 'none', 'DisplayName', 'Broadband range');
xlabel(ax2, 'Frequency [Hz]'); 
% ylabel(ax4, 'Power \muV^2 / Hz');
title(ax2, 'Average PSDs of successful trials')
grid(ax2, 'on');
xlim([0 350])

legend(legendUnq)
linkaxes([ax1 ax2], 'y')
fontsize(hFig, 16, 'points');

%% Plot Exemplary trials 
% Plot Raw signal, TF Map, Normalizer Output, and Cursor Position
ch = 1; % Channel to plot
trial= 8;  % Trial to plot
normOutputIdx = 1; % Select which Normalzier channel was used to control BCI task 
outField = sprintf('SGNL_NormalizerOutput_%d',normOutputIdx);
LCFreqs = [72.5 112.5];

% Baseline PSD
basePSD = data.baseline.RerefSignalPSD;

% Spectral parameters for TF map
winLength = round(0.25 * SR);              % 0.25 s window
winLength = 2^nextpow2(winLength);
wnd = hann(winLength, 'periodic');
noverlap = floor(0.8 * winLength);         % high temporal overlap
nfft = data.baseline.SpectralInfo.nfft;              % 1 Hz frequency resolution

% Smoothing and viusalization for TF map
fsmth = 10; % Frequency smoothing factor
tsmth = 10; % Time smoothing factor
f_cut = 400; % Cut-off frequency for visualization

% Convert Trial number to Run and Trial index
runTrials = cumsum(nTrials);
runIdx = find(trial <= runTrials, 1, 'first'); % Determine which run the trial belngs to
if runIdx == 1
    trialIdx = trial;
else
    trialIdx = trial-runTrials(runIdx-1);
end

BCITask = data.(fld_names{runIdx}).BCITask;
SR = data.(fld_names{runIdx}).SamplingRate;
trIdx = BCITask.startIdx(trialIdx):BCITask.endIdx(trialIdx);

% Get the trial signals
trialSignal = data.(fld_names{runIdx}).RerefSignal.Signal(trIdx,normOutputIdx);
normOutput = data.(fld_names{runIdx}).BCI2000ProcessingOutputs.(outField)(trIdx);
cursorPos = BCITask.CursorPosY(trIdx)./4095;

% Automatically detect LC bins
% LCField = 'Out2'; % Select which LCclassifier output to plot 
% LCFreqs = data.(fld_names{runIdx+1}).BCI2000ProcessingOutputs.classifierBins.(LCField).inputElementNumeric;

% Create figure
figure('Position', [200 100 1200 800])
tl = tiledlayout(3,1, ...
    'TileSpacing', 'compact', ...
    'Padding', 'compact');

t_ax = trIdx./SR; % Time axis

% Plot raw signal
ax1 = nexttile;
plot(t_ax, trialSignal, 'Color', [0.1 0.1 0.1], 'LineWidth', 1.5)
ylabel('Amplitude [µV]')
title(sprintf('Raw Signal, Ch %d', ch))
xlim([t_ax(1) t_ax(end)]);

% Time-frequency map
ax2 = nexttile;
[~, f, t, p] = spectrogram(trialSignal, wnd, noverlap, nfft, SR, 'power');
t_spec = t + t_ax(1);
% Normalize the power to baseline
p_norm = p./basePSD;
% Cut off the high freuqncies
mask = f>f_cut; f(mask)=[]; p_norm(mask,:)=[];

% Smooth TF map
p_norm = smoothdata(p_norm, 1, 'gaussian', fsmth); % smooth across frequency
p_norm = smoothdata(p_norm, 2, 'gaussian', tsmth); % smooth across time
% p_norm = imgaussfilt(p_norm, [0.6 0.3] );
imagesc(t_spec, f, p_norm); hold on
yline(LCFreqs, 'Color','m', 'LineWidth',2, 'LineStyle','--')
axis xy
cb = colorbar;
cb.Label.String = 'Normalized Power [n/a]';
ylabel('Frequency [Hz]')
title('Time-frequency map')
xlim([t_ax(1) t_ax(end)]);
ylim([0 180])
% clim(prctile(p_norm, [2.5 97.5], 'all'))
clim([0 10])
colormap(rgb.cm)
% colormap turbo

% 3) Normalized output and cursor position
ax1 = nexttile;
yyaxis left
plot(t_ax, normOutput, 'LineWidth', 1.5); hold on
yline(ax1, 0, 'k', 'HandleVisibility', 'off');
ylabel('Normalized Output')
% Set symmetric limits
currentLimits = ylim; 
maxLimit = max(abs(currentLimits));
ylim([-maxLimit, maxLimit]);

yyaxis right
plot(t_ax, cursorPos, 'LineWidth', 1.5)
ylabel('Cursor PosY')
xlabel('Time [s]')
title('Decoder Output and Cursor Position')
grid on
ylim([0 1])
xlim([t_ax(1) t_ax(end)])

% Link x axes
linkaxes([ax1 ax2 ax1], 'x')
sgtitle(sprintf('Trial %d', trial))

clear ch trial normOutputIdx outField LCField