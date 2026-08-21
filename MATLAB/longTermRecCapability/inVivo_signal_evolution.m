% INVIVO_SIGNAL_EVOLUTION - inVivo signal evolution
%
% Author: Frederik Lampert
% Institution: Mayo Clinic
% Year: 2026
%
% This script reproduces an in vivo analysis workflow used for evaluation
% of long-term neural recording performance and signal evolution.
%
% Data availability: source recordings are available in BCI2000 .dat format
% through OpenNeuro (doi:10.18112/openneuro.ds004624.v3.0.0).
%
% Before running:
%   1. Replace placeholder paths with local project/data locations.
%   2. Add the BIC_data class and required helper functions to the MATLAB path.
%   3. Review subject-specific selections, exclusions, and plotting settings.
%
% MATLAB version: R2023b

clear all; close all; clc
%% Configure project path and Add Paths
% Replace with the local repository or source-code folder.
projectPath = '/path/to/project/';
addpath(genpath(projectPath));
addpath(genpath('/path/to/BCI200mexfiles/'));

%% Specify constants and paths
SAVE_FIG = 1;
SAVE_STRUCT = false;
fig_pth = '/path/to/figure_output_directory/';
AMP_FACTORS = {'57.5dB','51.5dB','45.5dB','39.5dB'};    % set list of amplifaction factors for Cortec

%% Select directory with data and get filepaths
sub = 'sub-cXX';
data_folder = ['/path/to/data/' filesep sub];
recs = dir(fullfile(data_folder,'*.dat'));
[~, idx] = sort([recs.datenum]);   % ascending (oldest → newest)
recs = recs(idx);

% Filter out every N-th recording
N = 10;
recs = [recs(1:N:end); recs(end)];

switch sub
    case 'sub-c00' % Laika
        cnstrct = 'surface'; spec = 'NN';  % Custom difference pair for cortec Grids
        knownBadChnls = [];
        channelLeadMap = dictionary( ...
            'left', struct('indices', 1:16, 'labels',  "Ch" + (1:16)), ...
            'right',struct('indices', 17:32, 'labels',  "Ch" + (17:32)));
    case 'sub-c01' % Belka
        cnstrct = 'surface'; spec = 'NN';  % Custom difference pair for cortec Grids
        knownBadChnls = [21 22];
        channelLeadMap = dictionary( ...
            'frontal', struct('indices', 1:12, 'labels',  "Ch" + (1:12)), ...
            'temporal',struct('indices', 13:22, 'labels',  "Ch" + (13:22)), ...
            'occipital',struct('indices', 23:32, 'labels',  "Ch" + (23:32)));
        implantDate = datetime('March 20, 2023');
    case 'sub-c02'  % Billy
        cnstrct = 'depth'; spec = 'D-S-D-S';
        knownBadChnls = [2:7];
        channelLeadMap = dictionary( ...
            'LANT', struct('indices', 1:8, 'labels',  "Ch" + (1:8)), ...
            'LHPC',struct('indices', 9:16, 'labels',  "Ch" + (9:16)), ...
            'RANT',struct('indices', 17:24, 'labels',  "Ch" + (17:24)),...
            'RHPC', struct('indices', 25:32, 'labels',  "Ch" + (25:32)) );
        implantDate = datetime('May 13, 2024');
    case 'sub-c03'  % Strelka
        cnstrct = 'surface'; spec = 'NN';  % Custom difference pair for cortec Grids
        knownBadChnls = [11 21];
        channelLeadMap = dictionary( ...
            'frontal', struct('indices', 1:12, 'labels',  "Ch" + (1:12)), ...
            'temporal',struct('indices', 13:22, 'labels',  "Ch" + (13:22)), ...
            'occipital',struct('indices', 23:32, 'labels',  "Ch" + (23:32)));
        implantDate = datetime('Nov 04, 2024');
    case 'sub-c04' % Willie
        cnstrct = 'combined'; spec = 'G-G-D-S';
        knownBadChnls = [1:8];
        channelLeadMap = dictionary( ...
            'displaced', struct('indices', 1:8, 'labels',  "Ch" + (1:8)), ...
            'occipital',struct('indices', 9:16, 'labels',  "Ch" + (9:16)), ...
            'thalamus',struct('indices', 17:24, 'labels',  "Ch" + (17:24)),...
            'lobus-piriformis', struct('indices', 25:32, 'labels',  "Ch" + (25:32)) );
        implantDate = datetime('Jul 28, 2025');
    case 'sub-c05' % Mushka
        cnstrct = 'combined'; spec = 'G-G-G-D';
        knownBadChnls = [];
        channelLeadMap = dictionary( ...
            'frontal', struct('indices', 1:8, 'labels',  "Ch" + (1:8)), ...
            'temporal',struct('indices', 9:16, 'labels',  "Ch" + (9:16)), ...
            'occipital',struct('indices', 17:24, 'labels',  "Ch" + (17:24)),...
            'thalamus', struct('indices', 25:32, 'labels',  "Ch" + (25:32)) );
        implantDate = datetime('January 12, 2026');
end

%%  ----------------- Pre-processing parameters  -----------------

% ~~~~~~~~ Select re-referencing ~~~~~~~~~~~~~~~~~
reref = 'Native';
% reref = 'CAR';         % Common average
% reref = 'BP'; 
% reref = 'leadSpecific'; RR_chnls = [1 9 17 25]
% reref = 'leadCAR';   % Custom Common Average for cortec Grids

% ~~~~~~~~~~ Filtering Parameters ~~~~~~~~~~~~~~~~~
notch = false; % Specify whether notch filter should be used
stopWidth = 5;      
FILT = true;       % Band-pass filtering
High_cut = 0;      % Specify High Cut-off Frequency (in Hz) for Low-Pass filtration
Low_cut = 0.5;     % Specify Low Cut-off Frequency (in Hz) for High-Pass filtration

% % ~~~~~~ Spectral analysis parameters ~~~~~~~
wndw = 5; % Time window in seconds for spectral calculation
freqs = 0:1:500; % Frequencies for spectral calculation
time_window = 10;   % Length of the extracted signal segment
data = struct();

%%  ----------------- Process data  -----------------
for i = 1:length(recs)
    fpath = fullfile(recs(i).folder,recs(i).name);
    % Load data
    [signal, states, parameters, total_samples ] ...      
        = load_bcidat(fpath, '-calibrated');         % load using MEX package   
    signal = interpolate_lost_samples(signal,states,'linear');  % linear interpolation of lost samples
    
    % Extract filename
    fname =strsplit(recs(i).name,'.');                         % Extract the filename without extension
    fname = matlab.lang.makeValidName(fname{1}); % Generate a valied struct fieldname

    % Extract parameters
    FS = parameters.SamplingRate.NumericValue;              % Get sampling frequency (rate)
    ref_ch = parameters.ReferenceCh.NumericValue;           % Get reference channel
    recording_date = datetime(parameters.StorageTime.Value);
    if isempty(recording_date) % fallback if could not be read from BCI2000 parms
        recording_date = datetime(recs(i).date);
    end

    try 
        amp_f = AMP_FACTORS{parameters.AmplificationFactor.NumericValue + 1}; % Get the amplification factor
    catch
        amp_f = '';
    end

    % Number of samples in requested window
    win_samp = round(time_window * FS);
    % Time axis for plotting
    t_ax = (0:win_samp-1) / FS;

    % calculate PSD
    win_fxn=hann(wndw*FS); % windowing function
    noverlap=floor(length(win_fxn) * 0.5); %overlap between calculations
    [pxx,f] = pwelch(signal,win_fxn,noverlap,freqs,FS); % mean power spectrum across whole experiment

    %% Clip and filter the signal
    % ~~~~ Filter the data ~~~~~
    if FILT
        % ---- IIR high-pass from MNL lab -----
        % filt_sig = ieeg_butterhighpass(signal,Low_cut,FS); % High-pass filter data using 3-ord Butter
        % filt_sig = ieeg_butterlowpass(signal,High_cut,FS); % Low-pass filter data using 3-ord Butter
        % filt_sig = ieeg_butterpass(proc_signal, [Low_cut High_cut], FS, true); 

        % ---- IIR high-pass from CNEL lab -----
        filt_sig = fastButterFilt(signal, FS, 'high', Low_cut, 3, false);

        % ---- FIR high-pass filter using EEGLAB  -----
        % filt_sig = eegfilt(signal',FS,Low_cut,High_cut); % EEGLAB expects data as [channels x samples]
        % filt_sig = filt_sig';

        data.(fname).Preprocessing.HighCut = High_cut;     % Store filter parameters
        data.(fname).Preprocessing.LowCut = Low_cut;     % Store filter parameters
    end

    if notch        % If notch  filter is required
        filt_sig = ieeg_notch(filt_sig, FS, 60, [59 61]);
        data.(fname).Preprocessing.Notch = notch;  % Store if notch was used parameters
    end

    % ~~~~~ Clip the signal ~~~
    % Define preferred start range: between 1/4 and 3/4 of the recording
    lower_idx = round(total_samples / 4);
    upper_idx = round(3 * total_samples / 4);
    % Constrain upper bound so the full window stays inside the signal
    upper_idx = min(upper_idx, total_samples - win_samp + 1);
    % Safety check
    if upper_idx < lower_idx
        error('Recording is too short for the requested time window.');
    end
    % Draw random valid start index
    start_idx = randi([lower_idx, upper_idx]);
    % End index
    end_idx = start_idx + win_samp - 1;
    filt_sig = filt_sig(start_idx:end_idx, :);

    %% Store the data
    data.(fname).FilePath = fpath;    % Store file name
    data.(fname).FileName = fname;            % File name without extension
    data.(fname).Date = recording_date;
    data.(fname).Signal = filt_sig;    % Store filtered signal
    data.(fname).SamplingRate = FS;         % Store Sampling Rate 
    data.(fname).Reference = ref_ch;        % Store reference channel
    data.(fname).AmpF = amp_f;              % Store amplification factor
    data.(fname).PSD = pxx;
    data.(fname).SpectrafFrequencies = f;
end
clear i noverlap total_samples f freqs wndw AMP_FACTORS amp_f fldnames recording_date

%% Save results
% Save outputs to a repository-independent local results directory. pre-proccessed data
if SAVE_STRUCT == true
    save(fullfile(data_folder, [sub '.mat']), 'data','-v7.3');
else 
    disp('File will NOT be saved')
end

%% Plot results
% Plotting parameters below reproduce the analysis figures and may be adapted as needed. the data
if ~exist('data','var')
    load(fullfile(data_folder, [sub '.mat']));
end
N_ch = 1;
FS = 1000;
freqs = 0:1:500; % Hardcoded, can be changed to get this from the data Struct
recs = fieldnames(data);

recording_dates = [cellfun(@(r) data.(r).Date, recs)];
recording_dates = dateshift(recording_dates, 'start', 'day');
daysSinceImplant = days(recording_dates - implantDate);
calDur = between(implantDate, recording_dates);
recording_date_labels = string(recording_dates, 'dd-MMM-yyyy');
daysSinceImplant_labels = string(daysSinceImplant);

time_window = 10;   % Length of the extracted signal segmen
t_ax = linspace(0,time_window, time_window*FS);   % Creating time axis for plotting signal
N_sig = numel(recording_dates);
amp_multiplier = 50;

% ColorMap
off = 5;
cmap = flipud(plasma(N_sig+off)); 
cmap(1:off,:) = []; %get rid of the white color

% Iterate over each channel
for ch = 1:N_ch
    % Pre allocate the signal mat (t x N dates)
    ch_sig = nan(length(t_ax),N_sig);
    % Pre allocate the signal mat (t x N dates)
    ch_PSD = nan(length(freqs),N_sig);
    % aggregate data for each channel across all recordings
    for i = 1:length(recs)
        ch_sig(:,i) = data.(recs{i}).Signal(:,ch);
        ch_PSD(:,i) = data.(recs{i}).PSD(:,ch);
    end

    % Plot signal evolution 
    plotSignalOverTime(ch_sig, FS, amp_multiplier,daysSinceImplant_labels, cmap, sprintf('Channel %d',ch))
    if SAVE_FIG
        figname = [fig_pth filesep 'Ch' num2str(ch) '_signal_evolution' '_subset'];
        print(gcf,figname,'-dpng','-r300')
        print(gcf,figname,'-dsvg','-vector')
    end

    % Plot PSD evolution
    plotPSDOverTime(ch_PSD, freqs,daysSinceImplant_labels, cmap, sprintf('Channel %d',ch))
    if SAVE_FIG
        figname = [fig_pth filesep 'Ch' num2str(ch) '_PSD_evolution'];
        print(gcf,figname,'-dpng','-r300')
        print(gcf,figname,'-dsvg','-vector')
    end
end

%% 
function plotSignalOverTime(signals, FS, amp_multiplier,dates, colors, fig_title)
    nSigs = size(signals,2);
    N = size(signals,1);
    t_ax = linspace(0,N/FS,N);

    y_spacing = 20000;

    if mod(nSigs,2)==1
        shift_s=linspace(0-floor(nSigs/2),0+floor(nSigs/2),nSigs)*y_spacing;
    else
        shift_s=linspace(0-floor(nSigs/2),0+floor(nSigs/2)-1,nSigs)*y_spacing;
    end
    signals = signals*amp_multiplier;
    shifted_signal = signals-repmat(shift_s,N,1);
    
    f = figure('Position', [100 100 1500 800]);
    main_ax = axes('Parent', f);
    hold(main_ax, 'on')

    for i = 1:nSigs
            plt = plot(main_ax, t_ax,shifted_signal(:,i), ...
                'DisplayName', dates(i), ...
                'Color',colors(i,:)); 
            plt.DataTipTemplate.DataTipRows(end+1) = dataTipTextRow('Date',repelem(dates(i),N));
    end
    
    % Compute major and minor tick locations
    yticks_major = flip(shift_s) * (-1);
    % yticks_minor = sort([yticks_major - 50*amp_multiplier, yticks_major + 50*amp_multiplier]);
    ylimits = [min(shifted_signal(:))-y_spacing/2, max(shifted_signal(:))+y_spacing/2];

    % LEFT Y-axis: Channel names
    main_ax.YLim = ylimits;
    main_ax.YTick = yticks_major;
    main_ax.YTickLabel = flipud(dates);
    main_ax.YLabel.String = 'Days since Implant';

    % Colorbar
    colormap(main_ax, flipud(colors))
    cb = colorbar(main_ax);
    nTicks = 11; % Number of tick labels to display
    dateIdx = round(linspace(1, nSigs, nTicks)); % Corresponding dates
    cb.TickLabels = flipud(dates(dateIdx));
    % cb.Label.String = 'Recording date';
    cb.Label.String = 'Days since Implant';

    % X-axis
    main_ax.XLim = [min(t_ax), max(t_ax)];
    main_ax.XLabel.String = 'Time [s]';    
    
    if exist('fig_title', 'var')
        title(fig_title,'Interpreter','none');
    else
        title('Raw signal');
    end
    
    grid on;
    fontsize(16, "points")

    % -------- Draw scale bar (bottom right corner) ------
    % Define scale sizes
    scale_amp = 100 * amp_multiplier;
    scale_time = 1;  % seconds

    % Position: bottom-right, just below lowest trace
    x0 = main_ax.XLim(2) - scale_time;
    y0 = min(shifted_signal(:)) - 0.3 * y_spacing;

    % Vertical line (amplitude scale)
    line(main_ax, [x0, x0], [y0, y0 + scale_amp], 'Color', 'k', 'LineWidth', 2, 'HandleVisibility', 'off');
    text(main_ax, x0 + 0.05 * scale_time, y0 + scale_amp/2, '100 \muV', ...
        'FontSize', 10, 'VerticalAlignment', 'middle');

    % Horizontal line (time scale)
    line(main_ax, [x0, x0 + scale_time], [y0, y0], 'Color', 'k', 'LineWidth', 2, 'HandleVisibility', 'off');
    text(main_ax, x0 + scale_time/2, y0 - 0.05 * y_spacing, ...
        sprintf('%.1f s', scale_time), ...
        'FontSize', 10, 'HorizontalAlignment', 'center');
    % legend('Location', 'eastoutside');

    axis square
end

%%
function plotPSDOverTime(PSDs, freqs,dates, colors, fig_title)
% Plot
figure('Name', 'PSD Evolution', 'Position',[100 100 1500 800]);
hold on;
for i = 1:size(PSDs,2)
    if i==size(PSDs,2) % Plot the last line more thick
        plot(freqs, ...
            smoothdata(pow2db(PSDs(:,i)),'movmean',10), ...
            'LineWidth', 3, ...
            'Color', colors(i,:), ...
            'DisplayName', dates(i));
    else
        plot(freqs, ...
            smoothdata(pow2db(PSDs(:,i)),'movmean',10), ...
            'LineWidth', 1.2, ...
            'LineStyle','-',...
            'Color', colors(i,:), ...
            'DisplayName', dates(i));
    end
end
hold off;
grid on;
xlabel('Frequency (Hz)');
ylabel('Power [dB]');
if exist('fig_title', 'var')
    title(fig_title,'Interpreter','none');
else
    title('PSD Evolution');
end

% Colorbar
colormap(flipud(colors))
cb = colorbar();
nTicks = 11; % Number of tick labels to display (dafult se to 11)
dateIdx = round(linspace(1, numel(dates), nTicks)); % Corresponding dates
cb.TickLabels = flipud(dates(dateIdx));
% cb.Label.String = 'Recording date';
cb.Label.String = 'Days since Implant';

% legend('Location', 'eastoutside');
axis square
fontsize(14,'points')
end
