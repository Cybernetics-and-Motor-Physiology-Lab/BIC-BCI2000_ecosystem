clear all; close all; clc
addpath('path/to/BCI2000/MATLAB/mex');
addpath('path/to/mnl_ieegBasics/functions'); % https://github.com/MultimodalNeuroimagingLab/mnl_ieegBasics

%% ~~~~~~ Select paths ~~~~~~~~~~~~~
dtpth = 'pat/ro/data/';
sessions = {'shorted_channels_39dB_ref32001', 'shorted_channels_45dB_ref32001', 'shorted_channels_51dB_ref32001', 'shorted_channels_57dB_ref32001'};

fig_dir = '/path/to/figure_output_directory';    % specify path where to save the figures
save_dir = '/path/to/data_output_directory';
%% ~~~~~~~~~~ SET FLAGS constants & parameters~~~~~~~~~~
SAVE_SPECT = true;
SAVE_FIG = 0; %false;
dB_AXIS = false;

sp = filesep;
amps = {'57.5dB',  '51.5dB', '45.5dB',  '39.5dB'};      % Specify available amplification values
chnls = 1:32;                                           % chanels selected
% ----- Spectral parameters ------
FS = 1000;
freqs = 0:0.5:500;                                      % Specify frequency range for spectra
win_fxn='hann';                                         % windowing function
win_len = 2;                                           % Window length in seconds
noverlap_coeff = 0.9;                                         % Overlap cofficient
bands = dictionary( ...
    ["Delta", "Theta", "Alpha", "Beta", "LB", "HB", "Total"], ...
    {[1 4], [4 8], [8 15], [15 30], [1 80], [80 500], [2.5 500]}); % Specify frequency bands
clrs = {[0 0.4470 0.7410], [0.9290 0.6940 0.1250], [0.4660 0.6740 0.1880], [0.8500 0.3250 0.0980]}; % collors for plotting
transp = 0.1;

%% ~~~~~~~~~ Plot data ~~~~~~~~~~~~~~~~~
gain = 100;
BP_freqs = [1, 499]; % Band passiing only for visualization
for ses = 1:size(sessions,1)
    for r = 1:size(sessions,2)
        fl = get_paths([dtpth sp sessions{ses,r}],'.dat','first'); % select file from the session
        [signal, states, parameters, total_samples, file_samples ] ...
            = load_bcidat(fl, '-calibrated');  % load using MEX package
        fl_name = strsplit(fl,'/'); fl_name = strsplit(fl_name{end},'.'); fl_name = fl_name{1};
        % Band Pass signal (only for visualization)
        signal = ieeg_butterpass(signal, BP_freqs, parameters.SamplingRate.NumericValue, true); % High-pass filter data using 3-ord Butter
        % Plot the raw signal
        displayChannels(signal, parameters.SamplingRate.NumericValue, parameters.ReferenceCh.NumericValue, 1:10,gain, fl_name);
        if SAVE_FIG
            print(gcf,[fig_dir sp 'raw_signal_' fl_name '.svg'],'-dsvg','-r100','-vector');
            % saveas(gcf,[fig_dir sp 'raw_signal_' fl_name '.fig']);
            saveas(gcf,[fig_dir sp 'raw_signal_' fl_name '.png']);
        end
    end
end

%% ~~~~~~~~~ Calculate specrtra for recordings ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~`
sig_spectra = struct();         % create empty structy for spectra

for ses = 1:size(sessions,1)
    for r = 1:size(sessions,2)
        fl = get_paths([dtpth sp sessions{ses,r}],'.dat','first'); % select file from the session

        [signal, states, parameters, total_samples, file_samples ] ...
            = load_bcidat(fl, '-calibrated');  % load using MEX package

        packetLoss = calculatePacketLoss(states);               % get packet loss
        FS = parameters.SamplingRate.NumericValue;              % Get sampling frequency (rate)
        ses_name = strsplit(sessions{ses,r},'_');                 % Extract name of the medium from the name

        % Windowing setup
        switch lower(win_fxn)
            case 'hann'
                wndw_fxn = hann(win_len*FS, 'periodic');
            case 'hamming'
                wndw_fxn = hamming(win_len*FS, 'periodic');
            otherwise
                error('Unsupported window function: %s', win_fxn);
        end
        noverlap = floor(length(wndw_fxn) * noverlap_coeff);

        % Calculate the PSDs
        [mps,freqs_out] = pwelch(signal,wndw_fxn,noverlap,freqs,FS); % mean power spectrum across whole experiment
        binWidth = mean(diff(freqs_out));

        % Save data and info into struct
        sig_spectra.(ses_name{1}).(['A' ses_name{3}]).spectra = mps;
        sig_spectra.(ses_name{1}).(['A' ses_name{3}]).PacketLoss = packetLoss;
        sig_spectra.(ses_name{1}).(['A' ses_name{3}]).AmplFactor = amps{parameters.AmplificationFactor.NumericValue+1};  % Get amplification factor
        sig_spectra.(ses_name{1}).(['A' ses_name{3}]).ReferenceCh = parameters.ReferenceCh.NumericValue;           % Get refeence channel;
        sig_spectra.(ses_name{1}).(['A' ses_name{3}]).Frequencies = freqs_out;
        sig_spectra.(ses_name{1}).(['A' ses_name{3}]).FreqBands = bands;
        sig_spectra.(ses_name{1}).(['A' ses_name{3}]).FilePath = fl;

        % Caluclate and save RMS values
        sig_spectra.(ses_name{1}).(['A' ses_name{3}]).TotalRMS = sqrt(sum(mps * binWidth));
        bandNames = keys(bands);
        for b = 1:numel(bandNames)
            bandName = bandNames{b};
            bandFreq = bands{bandName};
            idx = freqs >= bandFreq(1) & freqs <= bandFreq(2);
            sig_spectra.(ses_name{1}).(['A' ses_name{3}]).([bandName 'RMS']) = sqrt(sum(mps(idx, :) * binWidth));
        end
        % Extract and save impedance
        if isfile(strrep(fl,'.dat','_info0.txt'))
            imp = readcell(strrep(fl,'.dat','_info0.txt'));                % load the impedance measuerment
            sig_spectra.(ses_name{1}).(['A' ses_name{3}]).impedance = imp;
        end
    end
end

%~~~~~~~~ Save the spectral data ~~~~~~~~~~~~~
if SAVE_SPECT
    save([save_dir sp 'noiseFloor_estimate.mat'], 'sig_spectra', '-v7.3');
end

%% Double axis PSD plot
f = figure('Position', [100 100 1000 800]);
sessions = fieldnames(sig_spectra);

t = tiledlayout(numel(sessions), 1, 'Padding', 'loose', 'TileSpacing', 'compact');
sgtitle('Shorted input recordings PSDs');
y1 =1e-3; y2= 1e4; % Y limits

for ses = 1:numel(sessions)
    runs = fieldnames(sig_spectra.(sessions{ses}));

    nexttile; % one row per session
    hold on;

    ax1 = gca;  % left axis (linear scale)
    ax2 = axes('Position', ax1.Position, ...
        'YAxisLocation', 'right', ...
        'XAxisLocation', 'bottom', ...
        'Color', 'none', ...
        'XColor', 'k', 'YColor', 'k');

    for r = 1:numel(runs)
        ampl_f = sig_spectra.(sessions{ses}).(runs{r}).AmplFactor;
        spectra = sig_spectra.(sessions{ses}).(runs{r}).spectra;
        totalRMS = sig_spectra.(sessions{ses}).(runs{r}).TotalRMS;
        % !!!!!! HARDCODED SECTION !!!!!!!
        spectra(:,20) = [];         % remove bad channel
        totalRMS(:,20) = [];        % remove bad channel
        % !!!!!! END of HARDCODED SECTION !!!!!!!

        % --- Plot on linear scale (left Y axis)
        axes(ax1); hold on;
        for ch = 1:size(spectra, 2)
            plt = plot(freqs, spectra(:,ch), ...
                'Color', [clrs{r} transp], 'HandleVisibility', 'off');
            plt.DataTipTemplate.DataTipRows(end+1) = ...
                dataTipTextRow('Ch', repelem(ch, length(freqs)));
        end

        mean_lin = mean(spectra, 2);
        plt = plot(freqs, smooth(mean_lin), ...
            'Color', clrs{r}, 'LineWidth', 2, 'DisplayName', ampl_f);
        plt.DataTipTemplate.DataTipRows(end+1) = ...
            dataTipTextRow('Ch', repelem({'Average'}, length(freqs)));

        ylabel(ax1, 'Power \muV^2/Hz');
        set(ax1, 'YScale', 'log');
        ylim(ax1, [y1 y2]);

        % --- Plot on dB scale (right Y axis)
        axes(ax2); hold on;
        mean_dB = pow2db(mean_lin);
        % plot(freqs, mean_dB, 'Color', clrs{r}, 'LineWidth', 1, 'HandleVisibility', 'off');

        % Use dummy lines to scale dB axis correctly
        ylim(ax2, pow2db([y1 y2]));
        ylabel(ax2, 'PSD (dB/Hz)');
    end

    % Sync x-axis between both axes
    axis(ax1, 'square');
    axis(ax2, 'square');
    xlim(ax1, [0 500]);
    xlabel(ax1, 'Frequency [Hz]');
    grid on
    linkaxes([ax1, ax2], 'x');
    fontsize(16, 'points');
    % legend(ax1, legendUnq(f), 'Location', 'northeast');
    legend(ax1,runs)
    title(ax1, sessions{ses});
end

%% ~~~~~~~~~~~~~~~~~~~ Plot Noise densities (not the same as PSD) - recording medium + amp. factors~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% In similar way to Ayyoubi et al.
% load([save_dir sp 'noiseFloor_estimate.mat']);
bandNames = {'Total'};
bandColors = [
    0.9 0.9 1.0;   % light blue
    1.0 0.9 0.9;   % light red
    % 0.85 1.0 0.85; % light green
    % 1.0 1.0 0.7;   % light yellow
    % 0.85 0.95 0.95 % light cyan
    % 0.95 0.8 1.0;  % light purple
    ];

f = figure('Position', [100 100 1000 800]);
y1 =10^(-1.5); y2= 10^(0.5); % Y limits
hold on;
sessions = fieldnames(sig_spectra);
for ses = 1:numel(sessions)
    runs = fieldnames(sig_spectra.(sessions{ses}));
    % t = tiledlayout(1,3); % horizontal'
    for r = 1:numel(runs)
        ampl_f = sig_spectra.(sessions{ses}).(runs{r}).AmplFactor;
        % nexttile
        % subplot(1,3,m)%m
        mps = sig_spectra.(sessions{ses}).(runs{r}).spectra;   % Extract Power Spectral Density
        freqs = sig_spectra.(sessions{ses}).(runs{r}).Frequencies;  % Extract frequencies from saved data
        binWidth = mean(diff(freqs)); %Extract

        mps(:,20) = []; % Remove BAD CHANNNEl !!!!! DELETE THIS IF USING FOR OTHER ANALYSIS!!!!!

        % Plot indiviudal channels
        for ch = 1:size(mps,2)
            plt = plot(freqs,sqrt(mps(:,ch)),'color',[clrs{r} transp],'HandleVisibility','off');      % 'DisplayName',ampl_f,
            plt.DataTipTemplate.DataTipRows(end+1) = dataTipTextRow('Ch',repelem(ch,length(freqs)));
        end

        % Plot average PSD
        avg_asd = smooth(mean(sqrt(mps), 2));
        plt = plot(freqs,avg_asd, 'color',clrs{r}, 'LineWidth', 2,'DisplayName',ampl_f); %movmean(mean(spectraNSD,1),3)
        plt.DataTipTemplate.DataTipRows(end+1) = dataTipTextRow('Ch', repelem({'Average'},length(freqs)));


        % Set scales, labels, etc
        yscale log; ylabel('Noise Spectral Density [\muV/\surdHz]');
        ylim([y1 y2])
        xlim([0 500]); xlabel('Frequency [Hz]');

        % Draw patches under highest line
        if r==1
            for b = 1:numel(bandNames)
                bandName = bandNames{b};
                bandFreq = bands{bandName};
                Xcoord = freqs(freqs > bandFreq(1) & freqs <= bandFreq(2));
                Ycoord = avg_asd(freqs > bandFreq(1) & freqs <= bandFreq(2));
                area(Xcoord, Ycoord, ...
                    'FaceColor',bandColors(b,:), ...
                    'FaceAlpha', 0.5, ...
                    'EdgeColor', 'none',...
                    'HandleVisibility','off');  % Optional: for legend
                uistack(plt,'top')
            end
        end

        % Caluclate RMS voltage values
        idx = freqs > 2.5;
        totalRMS = mean(sqrt(sum(mps(idx,:) * binWidth, 1))); % calculate values above the 2.5Hz to omit the DC
        for b = 1:numel(bandNames)
            bandName = bandNames{b};
            bandFreq = bands{bandName};
            idx = freqs >= bandFreq(1) & freqs <= bandFreq(2);
            bandRMS_per_ch = sqrt( sum(mps(idx,:) * binWidth, 1) );  % [µV]
            sig_spectra.(sessions{ses}).(runs{r}).([bandName 'noiseRMS']) = bandRMS_per_ch;
        end
        % --- Add text with RMS
        xText = gca().XLim(1) + r*(0.1*diff(gca().XLim));
        yText = 10^(log10(gca().YLim(2)) - r*0.05*diff(log10(gca().YLim)));
        txt = sprintf('Total RMS = %.4g', mean(totalRMS));

        text(xText, yText, txt, ...
            'FontSize', 12, ...
            'EdgeColor', clrs{r}, ...
            'Margin', 5);


        % Set legend
        % legend(legendUnq(f))
        legend(runs)

    end
    sgtitle('Shorted input PSDs');
    fontsize(16,'points');
    axis square
end

%% ~~~~~~~~~~~~~~~~~~~ Plot log-log spectra + Power-law quantification~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% load([save_dir sp 'saline_test_noisefloor.mat']);
fitRange = [10 499];
exclude = 60:60:FS/2;   % or 60:60:480 for 60 Hz mains
opts = struct('excludeHalfWidthHz',1,'useLogResidual',true);

clrs = {[0 0.4470 0.7410], [0.9290 0.6940 0.1250], [0.4660 0.6740 0.1880], [0.8500 0.3250 0.0980]}; % collors for plotting
transp = 0.1;
f = figure('Position', [100 100 1000 800]);
hold on;
sessions = fieldnames(sig_spectra);
for ses = 1:numel(sessions)
    runs = fieldnames(sig_spectra.(sessions{ses}));
    % t = tiledlayout(1,3); % horizontal'
    for r = 1:numel(runs)
        ampl_f = sig_spectra.(sessions{ses}).(runs{r}).AmplFactor;
        spectra = sig_spectra.(sessions{ses}).(runs{r}).spectra;         % Extract spectra
        % spectra(:,20) = []; % !!!!!remove bad channel
        % Plot average PSD
        p = plot(freqs,smooth(mean(spectra,2)), 'color',clrs{r}, 'LineWidth', 1,'DisplayName',ampl_f); %movmean(mean(spectra,1),3)
        hold on;

        % Fit
        out = fitPowerLawNoiseFloor(freqs, mean(spectra,2), fitRange, exclude, opts);
        plot(freqs, out.Pfit, 'color',[clrs{r} 0.8], 'LineWidth', 2, 'LineStyle','--');

        % Set scales, labels, etc
        set(gca, 'XScale', 'log', 'YScale', 'log')
        ylabel('Power \muV^2/Hz')
        xlim(fitRange); xlabel('Frequency [Hz]');

        ax = gca;

        % Use logarithmic spacing for stable placement
        xText = 10^(0.9 + r*0.15*diff(log10(ax.XLim)));
        yText = 10^(-1.5 - r*0.05*diff(log10(ax.YLim)));

        txt = {
            sprintf('A = %.3g', out.A)
            sprintf('\\chi = %.3f', out.chi)
            sprintf('C = %.3g', out.C)
            };

        text(xText, yText, txt, ...
            'FontSize', 10, ...
            'FontName', 'Helvetica', ...
            'VerticalAlignment', 'top', ...
            'BackgroundColor', 'w', ...
            'EdgeColor', clrs{r}, ...
            'Margin', 5);

        % Set legend
        legend(legendUnq(f))
    end
    sgtitle('Saline PSDs');
    fontsize(16,'points');
    axis square
end


if SAVE_FIG
    print([gcf,fig_dir fl_id{1} '_grid' num2str(grd) '.svg'],'-dsvg','-r100','-vector');
    saveas(gcf,[fig_dir filesep 'noise_floor_quantification.fig']);
    saveas(gcf,[fig_dir filesep 'noise_floor_quantification.png']);
end

