clear all; close all; clc;
%% ------------------  SELECT FILE ---------------------
[signal, states, parameters, total_samples, file_samples ] = LoadBCI2kData(); 
%% ----------------- Specify Parameters -----------------
%%%%%%%%% Specify recording construct %%%%%%%
cnstrct = 'depth';

%%%%%%%%% General parameters  %%%%%%%%%%%%
bad_chnls = [];    % Specify Bad channels
chnls = 1:size(signal,2);

%%%%%%%%% Spectral analysis parameters %%%%%%%
FS = parameters.SamplingRate.NumericValue;              % Get sampling frequency (rate)
freqs= 0:0.5:500;               % Frequencies for spectral analysis
alp = 0.1;                % Specify transparency for PSD Plots

%% ~~~~~~~~~ Plot raw Band-passed data ~~~~~~~~~~~~~~~~~
gain = 50;
BP_freqs = [1, 499];
% Band Pass signal
filt_signal = ieeg_butterpass(signal, BP_freqs, parameters.SamplingRate.NumericValue, true); % High-pass filter data using 3-ord Butter
% Plot the raw signal
displayChannels(filt_signal(5000:35000,:), parameters.SamplingRate.NumericValue, parameters.ReferenceCh.NumericValue, 2,gain, 'Raw');
% clear filt_signal
%% Re-reference data
% CAR
[CAR_signal, CAR_names] = cortec_reref(signal,cnstrct,bad_chnls,'CAR');
mean_signal = mean(signal(:,setdiff(chnls,bad_chnls)),2);

% Diference pairBipolar re-refernecing within the lead
[BP_signal, BP_names] = cortec_reref(signal,cnstrct,bad_chnls,'BP','NN');

% Diference pairBipolar re-refernecing within the lead
[leadCAR_signal, leadCAR_names] = cortec_reref(signal,cnstrct,bad_chnls,'leadCAR');

% Lead specific reference 'leadMono'
[leadMP_signal, leadMP_names] = cortec_reref(signal,cnstrct,bad_chnls,'leadSpecific',[1 14 25]);

%% Examine the re-referenced signal variance
fig = figure('Name', 'Rereferenced signal variance', ...
            'Position', [100 100 1500 800], 'Units','normalized');
mrksz = 10; % markersize

% Plot Raw signal
h(1) = plot(chnls, var(signal,0,1), 'color', [0 0.4470 0.7410 alp], ...
    'Marker','o', 'MarkerSize',mrksz, 'MarkerFaceColor','#0072BD'); hold on
% Plot CAR signal
h(2) = plot(chnls, var(CAR_signal,0,1), 'color', [0.8500 0.3250 0.0980 alp], ...
    'Marker','o', 'MarkerSize',mrksz, 'MarkerFaceColor','#D95319');
% Plot BP signal
h(3) = plot(cellfun(@(x) mean(str2double(split(x, '-'))), BP_names), ...
     var(BP_signal,0,1),'color', [0.9290 0.6940 0.1250 alp], ...
     'Marker','^', 'MarkerSize',mrksz,'MarkerFaceColor','#EDB120');
% Plot lead CAR signal
h(4) = plot(chnls, var(leadCAR_signal,0,1), 'color', ...
    [0.4940 0.1840 0.5560 alp],'Marker','o', 'MarkerSize',mrksz,'MarkerEdgeColor','#7E2F8E');
% Plot lead Monopolar signal
h(5) = plot(chnls, var(leadMP_signal,0,1), 'color', ...
    [0.3010 0.7450 0.9330 alp],'Marker','diamond', 'MarkerSize',mrksz, ...
    'MarkerFaceColor','#4DBEEE');

xlabel('Channel no.'); xlim([0 32]); xticks(0:1:32)
ylabel('Variance'); ylim([0 15000]);
yscale log
grid on

title('Effect of re-referencing on channel variance', 'Interpreter','none')
lgnd = legend({'Raw','CAR','DP','leadCAR', 'leadDP', 'leadSpeciffic'},'FontSize',14);
% fontsize(lgnd,14,'points')
fontsize(16,'points')

%% Calculate spectra for each re-referencing type
raw_spectra = NoiseFloorSpect(signal,FS,freqs);    % calculate mean power spectra for each channels
CAR_spectra = NoiseFloorSpect(CAR_signal,FS,freqs);    % calculate mean power spectra for each channels
BP_spectra = NoiseFloorSpect(BP_signal,FS,freqs);    % calculate mean power spectra for each channels
leadCAR_spectra = NoiseFloorSpect(leadCAR_signal,FS,freqs);    % calculate mean power spectra for each channels
leadSpec_spectra = NoiseFloorSpect(leadMP_signal,FS,freqs);    % calculate mean power spectra for each channels
sig_spectra = struct('Native', raw_spectra, 'CAR', CAR_spectra, 'BP', BP_spectra, 'leadCAR', leadCAR_spectra, ...
    'leadSpeciffic', leadSpec_spectra, 'Frequencies', freqs, 'BadChannels', bad_chnls);
%% Plot average spectra
fig = figure('Name', 'Signal PSD', ...
    'Position', [100 100 1500 800]);
% Get rid of bad channels
chnls(bad_chnls) = [];

% Plot Raw signal
h(1) = plot(freqs,pow2db(mean(raw_spectra(chnls,:),1)), 'LineWidth',2); hold on
% Plot CAR signal
h(2) = plot(freqs,pow2db(mean(CAR_spectra(chnls,:),1)), 'LineWidth',2);
% Plot BP signal
h(3) = plot(freqs,pow2db(mean(BP_spectra,1)), 'LineWidth',2);
% Plot lead CAR signal
h(4) = plot(freqs,pow2db(mean(leadCAR_spectra(chnls,:),1)), 'LineWidth',2);
% Plot lead specific signal
h(5) = plot(freqs,pow2db(mean(leadSpec_spectra(chnls,:),1)), 'LineWidth',2);

grid on
title('Noise floor')
ylabel('Power (dB/Hz)');
xlim([0 freqs(end)]);
xlabel('Frequency (Hz)');
lgnd = legend({'Raw','CAR','BP','leadCAR', 'leadSpeciffic'},'FontSize',14);
fontsize(20,'points');
   

% %% Examine the re-referenced signal variance
% fig = figure('Name', 'Rereferenced signal variance', ...
%             'Position', [100 100 1500 800], 'Units','normalized');
% 
% mrksz = 100; % markersize
% 
% switch cnstrct
%     case 'depth'
%         % Plot Raw signal
%         h(1) = scatter(chnls, var(signal,0,1), mrksz,'o', 'filled'); hold on
%         % Plot CAR signal
%         h(2) = scatter(chnls, var(CAR_signal,0,1), mrksz,'o', 'filled');
%         % Plot BP signal
%         h(3) = scatter(cellfun(@(x) mean(str2double(split(x, '-'))), DP_names), ...
%             var(DP_signal,0,1),mrksz,'^','filled');
%         % Plot lead CAR signal
%         h(4) = scatter(chnls, var(leadCAR_signal,0,1),mrksz,'o','LineWidth',2);
%         % Plot lead Monopolar signal
%         h(5) = scatter(chnls, var(leadMP_signal,0,1),mrksz,'diamond', 'filled','LineWidth',2);
%     case 'surface'
%         % Plot Raw signal
%         h(1) = scatter(chnls, var(signal,0,1), mrksz,'o', 'filled'); hold on
%         % Plot CAR signal
%         h(2) = scatter(chnls, var(CAR_signal,0,1), mrksz,'o', 'filled');
%         % Plot lead CAR signal
%         h(3) = scatter(chnls, var(leadCAR_signal,0,1),mrksz,'o','LineWidth',2);
%         % Plot lead Monopolar signal
%         h(4) = scatter(chnls, var(leadMP_signal,0,1),mrksz,'diamond', 'filled','LineWidth',2);
% end
% 
% xlabel('Channel no.'); xlim([0 32]); xticks(0:1:32)
% ylabel('Variance'); ylim([0 15000]);
% yscale log
% 
% grid on
% 
% title('Effect of re-referencing on channel variance', 'Interpreter','none')
% lgnd = legend({'Raw','CAR','DP','leadCAR', 'leadDP', 'leadSpeciffic'},'FontSize',14);
% fontsize(lgnd,14,'points')
% fontsize(16,'points')
