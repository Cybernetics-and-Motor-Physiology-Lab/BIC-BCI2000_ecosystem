clear all; close all; clc;
% -------------------------------------------------------------------------
%% --------------------- 0.1 Extract Chirp Data ---------------------------
% -------------------------------------------------------------------------
    % SELECT FILE 
    [signal, states, parameters, total_samples, file_path ] = LoadBCI2kData();
    [fileDirectory,fileName, fileExt] = fileparts(file_path); 
    signal = double(signal);
    amp_f = {'57.5dB','51.5dB','45.5dB','39.5dB'};    % set list of amplifaction factors
    AMP_F = amp_f{parameters.AmplificationFactor.NumericValue + 1}; % Get the amplification factor
    FS = parameters.SamplingRate.NumericValue;              % Get sampling frequency (rate)
    REF_CH = parameters.ReferenceCh.NumericValue;           % Get reference channel(s)
    %% Parameters
    CH = 1; %Selected channel for analysis
    T = 250; % Time period of the chirp signal
    F = linspace(1,501,T*FS);    % Frequency sweep
    dsFactor = 10;  % Down Sample factor; adjust based on signal length vs plot width
    Vpp_in = 2;     % Inserted chirp amplitude (in mV) from signal generator
    %% Specify paths
    save_dir = 'path/to/data_output_directory';
    fig_dir = 'path/to/figure_output_directory';
    ampStr = regexprep(AMP_F, '\.\d+', ''); % for file Naming
    cmap = 'parula';
    SAVE_FIG=0;
    %% Display RAW data
    gain = 10;
    % HP signal
    signalForPlot = ieeg_highpass(signal, FS);
    displayChannels(downsample(signalForPlot, dsFactor), FS/dsFactor, REF_CH, [CH REF_CH], gain, 'Raw signal'); hold on
    if isfield(states, 'StimulusCode')  % if stimulation was on
        plot(states.StimulusCode)
    end
    clear signalForPlot
    %% Extract chirp periods
     [marked_times, marked_ids, chirp_idx] = segment_signal(detrend(signal(:,CH),2,1:FS:length(signal)), [], FS, fileName);
    
    %% Extract the chirp periods and their enevelope
    chirp = cell(numel(unique(marked_ids)),1);
    marked_idxs = round(marked_times.*FS);  % transform time indices into samples
    hlp =1;
    % Extract signals
    for i = 1:2:length(marked_idxs)-1
        chirp_i = signal(marked_idxs(i):marked_idxs(i+1),:); % Extract signal 
        chirp{hlp} = detrend(chirp_i,2,1:FS:length(chirp_i)); % Detrend the slow discharge induced by the signal generator
        hlp = hlp+1;
    end
    
    %% Save manually annotated signals
    ch_name = sprintf('Ch%d', CH);
    save([save_dir filesep fileName '.mat'], ...
        "chirp", "ch_name", "marked_ids", "marked_times", "marked_idxs", "FS", "AMP_F", "REF_CH", "fileName", "T", "Vpp_in");

% -------------------------------------------------------------------------
%% --------------------- 0.2 Analyze Chirp Data ---------------------------
% -------------------------------------------------------------------------

    %% Load the signal
    load('/path/to/extractedData/XXX.mat');
    
    %% Calculate enevlopes and pad the signal
    chirpEnvelopes = cell(numel(chirp),2); % 1st column Upper Env 2nd column Lower ENv
    % Calculate envelopes
    for j=1:numel(chirp)
        [chirpEnvelopes{j,1},chirpEnvelopes{j,2}] = envelope(chirp{j},2500,'analytic'); % Calculate envelope
        % [upEnv, lowEnv]
    end
    
    % Pad the signals to equal length
    maxLen = max(cellfun(@length, chirp));
    chirp = cellfun(@(x) [x(:); nan(maxLen - length(x), 1)], chirp, 'UniformOutput', false);
    chirpEnvelopes = cellfun(@(x) [x(:); nan(maxLen - length(x), 1)], chirpEnvelopes, 'UniformOutput', false);
    t = linspace(0,maxLen/FS, maxLen);t = linspace(0,maxLen/FS, maxLen);
    
    %% --------- Visualize ------------
    % Each chirp repetition
    f1 = figure('Position', [100 100 1200 800], 'Name', ch_name);
    tiledlayout('flow')
    for i = 1:numel(chirp)
        nexttile
        % t = linspace(1,length(chirp{i})/FS, length(chirp{i}));
        plot(t,chirp{i}, Color='k'); hold on
        plot(t,smoothdata(chirpEnvelopes{i,1},'gaussian',30),'r-', 'LineWidth',2); %Upper enevelop 
        plot(t,smoothdata(chirpEnvelopes{i,2},'gaussian',30),'b--','LineWidth',2); % Lower envelope
        title(sprintf('n=%d',i))
        xlabel('Time [s]'); ylabel('Amplitude [µV]')
        xlim tight
    end
    sgtitle(['Amp. ', AMP_F, ' ', ch_name]);
    
    if SAVE_FIG
        print(gcf,[fig_dir filesep 'chirp_' ampStr '_' ch_name '.svg'],'-dsvg','-r100','-vector')
        saveas(gcf,[fig_dir filesep 'chirp_' ampStr '_' ch_name '.fig']);
    end
    
    %% TF map of each chirp repetition
    f1 = figure('Position', [100 100 1200 800], 'Name', ch_name);
    tiledlayout('flow')
    for i = 1:numel(chirp)
        nexttile
        % t = linspace(1,length(chirp{i})/FS, length(chirp{i}));
        spectrogram(rmmissing(chirp{i}),1024,948,1024,FS,'yaxis');
        % [~,fr,tr,pxx] = spectrogram(rmmissing(chirp{i}),1024,948,1024,FS,'yaxis');
        % imagesc(tr,fr,pow2db(pxx+eps));     
        % axis xy
        % xlabel('Time [s]'); ylabel("Frequency [Hz]")
        title(sprintf('n=%d',i))
        colormap(cmap)
        clim([-30 40])
        xlim tight
    end
    sgtitle(['Amp. ', AMP_F, ' ', ch_name]);
    
    if SAVE_FIG
        print(f1,[fig_dir filesep 'TFmap_' ampStr '_' ch_name '.svg'],'-dsvg','-r100','-vector');
        saveas(f1,[fig_dir filesep 'TFmap_' ampStr '_' ch_name '.fig']);
    end
    
    %% PSD of each chirp repetition
    f2 = figure('Position', [100 100 1200 800], 'Name', ch_name);
    tiledlayout('flow')
    for i = 1:numel(chirp)
        nexttile
        % Recalculate using high number of nfft points instead of specifying
        % frequencies -> this will induce zero padding and hopefully reduces
        % scalloping
        [ps,freqs] = pwelch(rmmissing(chirp{i}),1024,948,0:0.5:501,FS); % mean power spectrum across whole experiment
        plot(freqs,pow2db(ps), 'Color','k', 'LineWidth',1.5)
        title(sprintf('n=%d',i))
        ylabel('PSD (dB/Hz)')
        xlabel('Frequency (Hz)')
        fontsize(16,'points')
        xlim tight
    end
    sgtitle(['Amp. ', AMP_F, ' ', ch_name]);
    
    if SAVE_FIG
        kjm_printfig([fig_dir filesep 'PSD_' ampStr '_' ch_name '.svg'],[12 12]);
        saveas(gcf,[fig_dir filesep 'PSD_' ampStr '_' ch_name '.fig']);
    end
    
    %% Average enevelope
    fig = figure('Position', [100 100 1200 800], 'Name', ch_name);
    hold on;
    % for i=1:size(chirpEnvelopes,1)
    %     % plot(t,chirp{i}, 'Color',[0, 0, 0, 0.3]); % Plot all the signals
    %     plot(t,smoothdata(chirpEnvelopes{i,1},'gaussian',100),'Color',[1, 0, 0, 0.65], 'LineWidth',1.2); %Upper enevelope
    %     plot(t,smoothdata(chirpEnvelopes{i,2},'gaussian',100),'Color',[0, 0.4471, 0.4471, 0.65], 'LineWidth',1.2); %Upper enevelope)
    % end
    upMean = mean(cell2mat(chirpEnvelopes(:,1)'),2, 'omitmissing');
    lowMean = mean(cell2mat(chirpEnvelopes(:,2)'),2, 'omitmissing');
    plot(t,smoothdata(upMean,'gaussian',100),'Color',[1, 0, 0], 'LineWidth',3);  % Average upper envelope 
    plot(t,smoothdata(lowMean,'gaussian',100),'Color',[0, 0.4471, 0.7412], 'LineWidth',3);  % Average lower envelope 
    % Patch signal
    x_ds = downsample(t, dsFactor);
    fill([x_ds, fliplr(x_ds)], [downsample(lowMean, dsFactor)', fliplr(downsample(upMean, dsFactor)')], 'k', ...
         'EdgeColor', 'none', 'FaceAlpha', 0.7);  % full black fill
    xlabel('Time [s]'); ylabel('Amplitude [µV]')
    xlim tight
    legend({'Average upper envelope', 'Average lower envelope', 'Average signal'})
    if SAVE_FIG
        print(gcf,[fig_dir filesep 'AVGchirp_' ampStr '_' ch_name '.svg'],'-dsvg','-r100','-vector')
        saveas(gcf,[fig_dir filesep 'AVGchirp_' ampStr '_' ch_name '.fig']);
    end
    
% -------------------------------------------------------------------------
%% ---------------- 0.3 Visualize Transfer-Function Characteristic --------
% -------------------------------------------------------------------------
    %% Frequency attenuation
    % Load data from the sinewave emasurements
    data_pth = '.../data/transferFucntion/TF_response_Vpp.xlsx';
    dataTable = readtable(data_pth, 'Sheet', 'A57dB');
    myFreqs = dataTable.Frequency_Hz;
    myAmps = dataTable.Vpp_muV./10^3;
    
    % Extract data from my chirp meausrements
    smoothUpMean = smoothdata(upMean,'movmean',100);
    smoothLowMean = smoothdata(lowMean,'movmean',100);
    measured_Vpp =  abs(diff([smoothUpMean,smoothLowMean],1,2));        %
    measured_Vpp = smoothdata(measured_Vpp,'movmean',1000);
    measured_Vpp_dB = 20*log10(measured_Vpp./(Vpp_in*10^3)); 
    
    freqPoints = F; 
    dx = abs(length(freqPoints)  - length(measured_Vpp_dB)); % How much needs to be cutted
    % trimLeft = floor(dx / 2);
    % trimRight = ceil(dx / 2);
    
    if length(freqPoints) >= length(measured_Vpp_dB)    % If average extracted chirp is shorter than the sweep time
        % freqPoints = freqPoints(1 + trimLeft : end - trimRight);
        freqPoints = freqPoints(1:end-dx);
    else
        % measured_Vpp_dB =measured_Vpp_dB(1 + trimLeft : end - trimRight,:);
        measured_Vpp_dB =measured_Vpp_dB(1:end-dx);
    end
    
    centerIdx = length(freqPoints);
    
    fig = figure('Position', [100 100 1200 800], 'Name', ch_name);
    hold on;
    plot(freqPoints, measured_Vpp_dB(1:centerIdx),  '-k',  'LineWidth',1.5); 
    % plot(freqPoints, flipud(measured_Vpp_dB(centerIdx+1:end)),  '-k',  'LineWidth',1); 
    plot(myFreqs, 20*log10(myAmps/2),'--*', 'LineWidth',1,'MarkerSize',5,...
        'MarkerEdgeColor', [1.0000    0.4118    0.1608]);
    xlim([0.5, 1000])
    set(gca,'XScale','log');
    xlabel('Frequency (Hz)');
    ylabel('Gain (dB)');
    title('Transfer Function of BIC');
    legend({'Frequency Sweep: 1-501 Hz over 250s', 'Sine Wave measurements'});
    grid on;
    box off;
    fontsize(16,'points')
    
    if SAVE_FIG
        print(gcf,[fig_dir filesep 'ampAtten_' ampStr '_' ch_name '.svg'],'-dsvg','-r100','-vector')
        saveas(gcf,[fig_dir filesep 'ampAtten_' ampStr '_' ch_name '.fig']);
    end
    
    % Save mean attenuation 
    % save([save_dir filesep fileName '_normAtten.mat'], 'meanAtten', 'freqPoints');