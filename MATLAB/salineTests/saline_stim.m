clear all; close all; clc
addpath('/Users/lampert.frederik/Documents/MATLAB/mex');
addpath('/Users/lampert.frederik/Documents/MATLAB/mnl_ieegBasics/functions');
addpath('/Users/lampert.frederik/Documents/MyCodes/MATLAB/functions');

%% ~~~~~~ select paths and load the data ~~~~~~~~~~~~~
fig_dir = '/Users/lampert.frederik/Documents/Mayo/CorTec/Figures/saline-test/implant04/stim/ses03';    % specify path where to save the figures
ddr = '/Users/lampert.frederik/Documents/Mayo/CorTec/data/other_recordings/Saline-tests'; % select data directory to open                          

%% ---------- Select file ----------------------
[fl, fl_pth] = uigetfile([ddr '/*.dat'], 'Select a File');                 % select file
fl_id = strsplit(fl, '.');                                                 % get file name and extension
[signal, states, parameters, total_samples, file_samples ] ...      
    = load_bcidat([fl_pth, fl], '-calibrated');         % load using MEX package 
data = interpolate_lost_samples(signal, states);        % interpolate lost values
%% ~~~~~~~~~~ Specify constants & Get parameters ~~~~~~~~~~~~~~~~
FS = parameters.SamplingRate.NumericValue;              % Get sampling frequency (rate)
T_ax = linspace(0,size(signal,1)/FS, size(signal,1));   % Creating time axis for plotting signal
amps = {'57.5dB',  '51.5dB', '45.5dB',  '39.5dB'};      % Specify available amplification values
ref_ch = parameters.ReferenceCh.NumericValue;           % Get refeence channel
sp = filesep;

%##############
SAVE_FIG = false;
%##############

chnls = [1:32];           % chanels selected

%% Plot the raw signal
displayChannels(data, FS, ref_ch, chnls, fl_id{1});
if SAVE_FIG
    % kjm_printfig([fig_dir fl_id{1} '_grid' num2str(grd)],[12 10]);
    saveas(gcf,[fig_dir sp fl_id{1} '_raw.fig']);
    saveas(gcf,[fig_dir sp fl_id{1} '_raw.png']);
end
%% ~~~~~~~~~~ Filter the signal ~~~~~~~~~~~~
% data = ieeg_notch(data, FS, 60, 2, true);            % Notch filter the data
% data = ieeg_butterhighpass(data,[0.05, 0.5],FS);      % High-pass filter data using 3-ord Butter

%% ~~~~~~~~~ Plot the chirp signal in time domain ~~~~~~~~~~~~~~~~~
chirp_start = find(diff(double(states.StimulusCode))==1);
chirp_stop = find(diff(double(states.StimulusCode))==-1);

grd = 1;
rep = 3;                 % select which repetition to plot

figure('Position', [100 100 1500 800]);
lgnd_lbls = {'Raw signal', 'Implant stimulation'};     % Labels for the legend
% clrs = ['r','y', 'b', 'c', 'm', 'g'];
clrs = ['r','r', 'r', 'r', 'r', 'r'];

for ch = 1:length(chnls)
    subplot(length(chnls),1,ch);            % create subplot
    plot(T_ax, data(:,chnls(ch)));        % Plot the raw signal for specific channel
    % xlim([0 T_ax(end)]);
    xlim([T_ax(chirp_start(rep)) T_ax(chirp_stop(rep))]);
    hold on
    title(sprintf('Channel %d',chnls(ch)), 'Interpreter', 'none');
    xlabel('Time [s]'); ylabel('U [uV]');
end

sgtitle({['Recording ', fl_id{1}, 'Amp.f ' amps{parameters.AmplificationFactor.NumericValue+1}],...
    ['  Ref.ch ', num2str(parameters.ReferenceCh.NumericValue)]}, ...
    'Interpreter', 'none') 
legend(lgnd_lbls);

if SAVE_FIG
    % kjm_printfig([fig_dir fl_id{1} '_grid' num2str(grd)],[12 10]);
    saveas(gcf,[fig_dir fl_id{1} 'ch2-20' num2str(grd) '.fig']);
    saveas(gcf,[fig_dir fl_id{1} 'ch2-20' num2str(grd) '.png']);
end

%% ~~~~~~~~~ Internal stimulation signal ~~~~~~~~~~~~~~~
stim_idx = find(diff(double(states.StimulusCode)));
% grd = 1;

stimulus_codes = unique(states.StimulusCode); stimulus_codes(1) = [];

% clrs = ['r','y', 'b', 'c', 'm', 'g'];
clrs = ['r','g', 'b'];

fig_handles = zeros(length(stimulus_codes),1);

for stm = 1:length(stimulus_codes)
    stim_period = signal(states.StimulusCode==stimulus_codes(stm),:);       % extract data nly for stimulation periods
    source_ch = parameters.StimulationTriggers.NumericValue(3, stimulus_codes(stm));
    dest_ch = parameters.StimulationTriggers.NumericValue(4, stimulus_codes(stm));
    stim_dur = parameters.StimulusDuration.NumericValue; % extract stimulus duration

    fig_handles(stm) = figure('Position',get(0, 'Screensize'));
    tlo = tiledlayout(ceil(length(chnls)/4),4,'TileSpacing','Tight');
    for ch = 1:length(chnls)
        % subplot(ceil(length(chnls)/4),4,ch);            % create subplot
        ax = nexttile;
        plot(stim_period(:,chnls(ch)));        % Plot the raw signal for specific channel
        xline(stim_dur*FS);xline(stim_dur*2*FS)
        xlim([0 length(stim_period)]);
        if ch == source_ch || ch == dest_ch
            ax.XColor = [1 0 0];
            ax.YColor = [1 0 0]; 
            title(sprintf('Channel %d',chnls(ch)), 'Interpreter', 'none', 'Color','r');
        else
            title(sprintf('Channel %d',chnls(ch)), 'Interpreter', 'none');
        end
        xlabel('Samples [-]'); ylabel('U [uV]');
    end

    sgtitle({['Recording ', fl_id{1}, ' Amp.f ' amps{parameters.AmplificationFactor.NumericValue+1}],...
    ['Stim. channels: ', num2str(source_ch), ...
    '-', num2str(dest_ch)...
    '  Ref.ch ',num2str(parameters.ReferenceCh.NumericValue) ]}, ...
    'Interpreter', 'none')


    if SAVE_FIG
        % kjm_printfig([fig_dir fl_id{1} '_grid' num2str(grd)],[14 12]);
        % saveas(gcf,[fig_dir fl_id{1} '_grid' num2str(grd).fig']);
        % saveas(gcf,[fig_dir sp fl_id{1} '_ch' num2str(source_ch) '-' num2str(dest_ch) '.fig']);
        saveas(gcf,[fig_dir sp fl_id{1} '_ch' num2str(source_ch) '-' num2str(dest_ch) '.png']);
    end
end


%% ------ Calulcate and plot the spectra ------
% ~~~~~ SPECIFY PARAMETERS FOR CALCULATION OF SPECTRA ~~~~~~~
win_fxn = hann(1*FS); % windowing function
noverlap = floor(length(win_fxn)/2); %overlap between calculations
freqs = 1:500;       % frequencies included in spectra calculation
% ~~~~~ Calculate spectra ~~~~~~~~
[mps,spect_f] = pwelch(data,win_fxn,noverlap,freqs,FS);  % mean power spectrum across whole experiment

for ch = 1:length(chnls)
    figure('Position', [100 100 1500 800])
    tiledlayout(2,1)
    nexttile
    semilogy(mps(:,ch),'r')
    xlabel('Frequency [Hz]'); ylabel('Power \muV^2/Hz')
    nexttile
    spectrogram(data(:,chnls(ch)),1*FS, noverlap, freqs, FS, 'yaxis');
    sgtitle(['Channel: ', num2str(chnls(ch)), ' ', fl_id{1}, ' Amp.f ',...
        amps{parameters.AmplificationFactor.NumericValue+1}], 'Interpreter', 'none')
    fontsize(16,'points');
    
    if SAVE_FIG
        % kjm_printfig([fig_dir fl_id{1} '_ch' num2str(chnls(ch))], [10,8]);
        % saveas(gcf,[fig_dir sp fl_id{1} '_ch' num2str(chnls(ch)) '.fig']);
        saveas(gcf,[fig_dir sp fl_id{1} '_ch' num2str(chnls(ch)) '.png']);
    end
end
