clear all; close all; clc
%% Add paths if necesscary
addpath('path/to/mnl_ieegBasics/functions'); % https://github.com/MultimodalNeuroimagingLab/mnl_ieegBasics
addpath(genpath('path/to/eeglab/functions')); % https://github.com/sccn/eeglab/tree/develop/functions
addpath(genpath('path/to/colormaps'));
bci2000path -AddAll

%% Load colormaps
rwb = load('rwb_colormap.mat');
rgb = load('RYGyCB.mat');

%% ------------------  SELECT FILE ---------------------
[signal, states, parameters, total_samples, file_path ] = LoadBCI2kData();
[~, fname, ~] = fileparts(file_path);
signal = double(signal);
raw_signal = signal; % store raw signal in case of reloading

if isfile(replace(file_path,'.dat', '_info0.txt'))
   impedanceValues = extractImpedanceValues(replace(file_path,'.dat', '_info0.txt'));
end

% Transmited channels for online-pocessig (account for 0-based indexing)
transChannels = cellfun(@(s) str2double(regexp(s, '\d+$', 'match', 'once')), parameters.TransmitChList.Value)+1;

%% Fallback if states were not dumped into SGNL fields
    % Calculate BCI2000 Output (Normalizer, etc)
    sigSrcFilters = 'TransmissionFilter';    %strjoin(parameters.SignalSourceFilterChain.Value(hlp,1), '|');
    % clear hlp;
    sigProcFilters = strjoin(parameters.SignalProcessingFilterChain.Value(1:4,1), '|'); % No control state filter
    filterChain = [sigSrcFilters '|' sigProcFilters];
    expr_output = bci2000chain(file_path, filterChain);
    % For complete output
    sigProcFilters = strjoin(parameters.SignalProcessingFilterChain.Value(1:end-1,1), '|'); % No control state filter
    filterChain = [sigSrcFilters '|' sigProcFilters];
    norm_output = bci2000chain(file_path, filterChain);

  % Housekeeping 
  clear expr_output sigProcFilters

%% --------- Select Subject ----------
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
    case 'c05'
        cnstrct = 'combined'; spec = 'G-G-G-D';
        knownBadChnls = [];
end

% ~~~~~ Recording parameters ~~~~~~~~~~~
AMP_FACTORS = {'57.5dB','51.5dB','45.5dB','39.5dB'};    % set list of amplifaction factors
amp_f = AMP_FACTORS{parameters.AmplificationFactor.NumericValue + 1}; % Get the amplification factor
FS = parameters.SamplingRate.NumericValue;              % Get sampling frequency (rate)
REF_CH = parameters.ReferenceCh.NumericValue;           % Get reference channel(s)
interpMethod = 'linear';

%% ----------------- SPECIFY Filtering and spectral Parameters, Specify Constants  -----------------
% ~~~~~~~~~~ Filtering Parameters ~~~~~~~~~~~~~~~~~`
% Select data filtration - if no filtration just comment out
BPf = [0.5 150];  % Specify Band Pass cut off frequency in Hz
BP_filt =true; % Specify whether band-pass filter should be used 
notch = false; % Specify whether notch filter should be used 

% ~~~~~~ Spectral analysis parameters ~~~~~~~
% Specify frequency range for PSD calculation
freqVec = 0.5:0.5:150;
% Windowing setup
wndw_fxn = hann(2^nextpow2(0.5*FS), 'periodic');        % Hann
% wndw_fxn = hamming(2^nextpow2(1*FS), 'periodic');   % Hamming
noverlap = floor(length(wndw_fxn) * 0.85);
nfft = length(wndw_fxn);
dB_AXIS = true;             % Choose y-axis labels
alpha = 0.8;                % Specify transparency for PSD Plots

% ~~~~~~ Statistical analysis parameters ~~~~~~~
pthresh=.05; % Specify p - treshold value

%% Re-reference data
[ref_sig, outputLabels, info] = applyBCI2kSpatialFilter(signal(:,transChannels), parameters, parameters.TransmitChList.Value); 

%% Mask the stim periods
% Specify time constant for for signal removal after a stimulation
stimArtifactT = 250; % in ms
% Mask stim periods
stim_onsets = find(diff(single(states.ImplantStimulation))==1);
requested_stim=find(diff(single(states.RequestedStimulation))==1);
mask = false(size(total_samples));

for j=1:length(stim_onsets)
    if stim_onsets(j)+stimArtifactT > size(signal,1)
        mask(stim_onsets(j)-10:size(signal,1)) = true;
    else
        mask(stim_onsets(j)-10:stim_onsets(j)+stimArtifactT-10) = true; % ten to get rid of the falling edge
    end
end

%Mask signal (get rid of stimulation artefact)
ref_sig(mask,:) = NaN;

%% Filter the signal
filt_sig = chunk_filter(ref_sig, FS, 'HP', BPf(1), 'LP', BPf(end), 'notch', [], 'order', 3);

%% Extract Spectral features
filt_sig = fillmissing(filt_sig, 'linear');     % Fill Nans

% Calculate spectograms
nCol = fix((total_samples-noverlap)/(length(wndw_fxn) - noverlap));                   % Get number of columns
spctgrms = complex(nan(length(freqVec),nCol,size(filt_sig,2)));     % Initiate empty matrix for spectograms
ps = zeros(length(freqVec),nCol,size(filt_sig,2));
for i=1:size(filt_sig,2)
    [spctgrms(:,:,i), f, t, ps(:,:,i)] = spectrogram(filt_sig(:,i), wndw_fxn, noverlap, freqVec, FS, "power");
end
Pdb = pow2db(ps + eps);

% ~~~~ Select power Bands based on the CL algorithm ~~~
% SPECTRAL ESTIMATION
  classifierBins = parseBCI2000Classifier(parameters.Classifier);
% PIB using IIR and Hilbert
    % classifierBins.(sprintf('Out%d',1)).binSpan = ...
    %     [parameters.HighPassCorner.NumericValue ... 
    %     parameters.LowPassCorner.NumericValue]; 

% Calculate envelope and filter the signal
HilbertPower = zeros(size(filt_sig));
HibertPhase = zeros(size(filt_sig));
for outCh = 1:size(filt_sig,2)
    freqBand = classifierBins.(sprintf('Out%d',outCh)).binSpan;
    [HilbertPower(:,outCh), HibertPhase(:,outCh)] = chunk_getHilbert(filt_sig(:,outCh), freqBand, FS, 'power');
end
% Control expressions
controlExpression = parameters.ControlExpression.Value{1, 1};
tokens = regexp(controlExpression, '(?:>=|<=|==|~=|>|<)\s*([-+]?\d*\.?\d+)', 'tokens');
thresholds = cellfun(@(t) str2double(t{1}), tokens);

clear freqBand tokens controlExpression
%% Segment task blocks
if isfield(states, 'StimulusCode')
    sc = single(states.StimulusCode);      % extracted from BCI2000 Stimulus Code
    [tc_dur, tc] = runlength(sc, length(sc));
    task_onsets = [1; cumsum(tc_dur(1:end-1)) + 1];
    taskTrials = find(tc == 1);
end

%% Prepre for visualization

% Smooth spectograms
Pdb_smooth = movmean(Pdb, 8, 1);   % smooth across frequency
Pdb_smooth = movmean(Pdb_smooth, 8, 2); % smooth across time

% Mask stim periods

HiblertP_smooth = LPFilter(HilbertPower, parameters.LPTimeConstant.NumericValue);
filt_sig(mask,:) = NaN;

% Remove frequency-wise background
% baseline = median(Pdb, 2);
% Pdb_smooth = Pdb_smooth - baseline;

%% Visualize whole signal
t_ax = linspace(0,total_samples/FS,total_samples);
ch = 1; % Specify which channels to plot
ds_factor = 4;

figure('Name','ClosedLoop', 'Position', [300 100 1200 800])
tiledlayout(4,1)

% Raw signal
nexttile

% Downsample signal
% sig_ds = smoothdata(resample(filt_sig(:,ch), 1, ds_factor),'sgolay',20);
% % Adjust time axis
% t_ds = t_ax(1:ds_factor:end);
% plot(t_ds,sig_ds, 'Color','#2b2b2b', 'LineWidth',1); hold on

% Original signal
plot(t_ax,filt_sig(:,ch), 'Color','#2b2b2b','LineWidth',1); hold on

xline(stim_onsets/FS, 'Color', [0.9294 0.6941 0.1255 alpha], 'LineWidth', 2)
% ylim(prctile(filt_sig(:,ch), [0.1 99.9]))
ylim([-200 200])
ylabel('Amplitude [µV]'); xlabel('Time [s]');
legend({'Filtered Signal', 'Stimulation'}, 'Location', 'northeastoutside')
xlim tight
title('Signal')

% Hilbert power
nexttile
plot(t_ax,HiblertP_smooth(:,ch), 'Color',[0.5, 0.5, 0.5], 'LineWidth', 1.5); hold on
xline(stim_onsets/FS, 'Color', [0.9294 0.6941 0.1255 alpha], 'LineWidth', 2)
ylabel('Power [µV]'); xlabel('Time [s]');
legend({'Hilbert Power', 'Stimulation'},'Location', 'northeastoutside')
xlim tight
title(sprintf('Hilbert Power within %d-%d Hz', classifierBins.(sprintf('Out%d',ch)).binSpan(1), classifierBins.(sprintf('Out%d',ch)).binSpan(2)))

% spectogram
nexttile
imagesc(t, f, Pdb_smooth(:,:,ch)); axis xy; hold on
xline(stim_onsets/FS, 'Color', [0.9294 0.6941 0.1255 alpha], 'LineWidth', 2)
colormap('turbo')
clim(prctile(Pdb_smooth(:), [2.5 99.5]));
cb = colorbar;
cb.Label.String = 'Power (dB/Hz)';
yline(classifierBins.(sprintf('Out%d',ch)).binSpan, 'Color',[1 0.0745 0.6510], 'LineWidth',1.5)
xlabel('Time (s)')
ylabel('Frequency (Hz)')


% Normalizer
nexttile
plot(t_ax, states.(sprintf('SGNL_NormalizerOutput_%d', ch-1)), 'Color', 'k','LineWidth', 2); hold on
% norm_offset = normalize(states.(sprintf('NormalizerOffset%d', ch)),'medianiqr');
% norm_gain = normalize(states.(sprintf('NormalizerGain%d', ch)),'medianiqr');
% plot(t_ax, norm_offset, 'Color', 'b','LineWidth', 1.2, 'Color',[0.9020 0.6824 0.6824])
% plot(t_ax, norm_gain, 'Color', 'g','LineWidth', 1.2, 'Color',[0.7882 0.9216 0.7529])
yline(thresholds(ch),'LineStyle','--', 'Color','r', 'LineWidth', 1.5)
xline(stim_onsets/FS, 'Color', [0.9294 0.6941 0.1255 alpha], 'LineWidth',2)
ylim(gca,[-5 5]);
ylabel('Normalizer Output [-]'); xlabel('Time [s]');
% legend({'Normalizer Output', 'Normalizer Offset', 'Normalizer Gain', 'Threshold','Stimulation'},'Location', 'northeastoutside')
legend({'Normalizer Output', 'Threshold','Stimulation'},'Location', 'northeastoutside')
xlim tight

ax = findall(gcf, 'Type', 'axes');   % get all axes in the figure
linkaxes(ax, 'x');                  % link x-axis (most common for tim

sgtitle(sprintf('%s: Output Ch %d',fname, ch),'Interpreter', 'none')
fontsize(16,'points');

clear t_ds %sig_ds norm_offset norm_gain
%% Visualize centered around the stim
t_w = 5; % time window for visualization
ch = 1; % Specify which channels to plot
t_ax = linspace(0,total_samples/FS,total_samples);
stimTrial = 4;

figure('Name','ClosedLoop', 'Position', [300 100 1200 800])
tiledlayout(3,1)

% Get t_indies centered around stim
t_idx = [stim_onsets-(t_w/2)*FS, stim_onsets+(t_w/2)*FS-1];
% Take care about theborderline values
t_idx(t_idx<1) =1;
t_idx(t_idx>total_samples) =total_samples;
currentWindow = t_idx(stimTrial,1):t_idx(stimTrial,2);
currentSpectralWindow = find(t<(t_idx(stimTrial,1)/FS),1,'last'):find(t>(t_idx(stimTrial,2)/FS),1,'first');

% Raw signal and Hilbert power
nexttile
plot(t_ax(currentWindow),filt_sig(currentWindow,ch), 'Color','#2b2b2b','LineWidth',1); hold on
plot(t_ax(currentWindow),HilbertPower(currentWindow,ch), 'Color',[0.65, 0.65, 0.65], 'LineWidth', 2)
xline(requested_stim(stimTrial)/FS, 'Color', [0.9294 0.6941 0.1255 alpha], 'LineWidth', 2, 'DisplayName','Requested Stimulation')
xline(stim_onsets(stimTrial)/FS, 'Color', [1 0 0 alpha], 'LineWidth', 2, 'DisplayName','Stimulation')
% ylim(prctile(filt_sig(currentWindow,ch), [2. 97.5]))
ylabel('Amplitude [µV]'); xlabel('Time [s]');
legend({'Signal', 'Hilbert Power', 'RequestedStimulation', 'Stimulation'})
axis tight

% Spectogram
nexttile
imagesc(t(currentSpectralWindow), f, Pdb_smooth(:,currentSpectralWindow,ch)); axis xy; hold on
colormap('turbo')
clim(prctile(Pdb_smooth(:), [5 99.5]));
cb = colorbar;
cb.Label.String = 'Power (dB/Hz)';
yline(classifierBins.(sprintf('Out%d',ch)).binSpan, 'Color',[1 0 0])
xlabel('Time (s)')
ylabel('Frequency (Hz)')
title(sprintf('Ch %d', ch))

% Normalizer
nexttile
chNormOutput = states.(sprintf('SGNL_NormalizerOutput_%d', ch-1));
plot(t_ax(currentWindow), chNormOutput(currentWindow), 'Color', 'k','LineWidth', 2); hold on
yline(thresholds(ch))
xline(requested_stim(stimTrial)/FS, 'Color', [0.9294 0.6941 0.1255 alpha], 'LineWidth', 2, 'DisplayName','Requested Stimulation')
xline(stim_onsets(stimTrial)/FS, 'Color', [1 0 0 alpha], 'LineWidth', 2)
ylabel('Normalizer Output [-]'); xlabel('Time [s]');
legend({'Normalizer Output', 'Threshold','RequestedStimulation', 'Stimulation'})
axis tight

clear norm_output
