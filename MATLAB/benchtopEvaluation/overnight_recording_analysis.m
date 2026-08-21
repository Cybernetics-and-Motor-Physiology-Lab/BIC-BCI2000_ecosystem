clear all; close all; clc;
%% Data Path
file_path = '/path/to/LongRecS003R01.dat';

%% Load Data
% Load just states
[~, states, parameters, total_samples] = load_bcidat(file_path, '-calibrated');

%% Convert states
% U = (double(states.FnirsiUSB_Voltage - 2^31)) .* 1e-5;
% I = (double(states.FnirsiUSB_Current - 2^31)) .* 1e-5;
P = (double(states.FnirsiUSB_Power   - 2^31)) .* 1e-5;
temperature = double(states.ImplantTemperature)./100;
lostSamples = double(states.ImplantLostSample);
% lostSamplePercentage = double(states.LostSampleBlockPercentage);

%% Downsample the data
fs = parameters.SamplingRate.NumericValue;
t_ax = 0:total_samples-1;
dsF = 100;
P = resample(P,1,dsF);
temperature = resample(temperature,1,dsF);
% lostSamplePercentage = downsample(lostSamplePercentage,dsF);
t_ax = t_ax(1:dsF:end)./fs;

%% Prepare for plotting
% Power
np_P = 500;  % adjust as needed
[upperP, ~] = envelope(P, np_P, 'peak');
upperP = smoothdata(upperP, 'SamplePoints', t_ax);

% Tmperature
tempSmooth = smoothdata(temperature, 'movmedian', 600);

%% PacketLoss 
PL_twindow = 60; % In seconds

% np_PL = 500;  % adjust as needed
% [upperLoss, ~] = envelope(lostSamplePercentage, np_PL, 'peak');
% upperLossSmooth = smoothdata(upperLoss, 'SamplePoints', t_ax);
avgPL = calculatePacketLoss(states); 

% Window length in samples
winSamp = round(PL_twindow * fs);
% Number of complete windows
nWins = floor(length(lostSamples) / winSamp);
% Trim excess samples so reshape works cleanly
trimLen = nWins * winSamp;
lostSamples_trim = lostSamples(1:trimLen);
% Reshape into windows
lostSamples_win = reshape(lostSamples_trim, winSamp, nWins);
% Percentage of lost samples per window
lostSamplePercentage = ...
    sum(lostSamples_win, 1) ./ winSamp * 100;
t_loss = ((0:nWins-1) * PL_twindow); % time axis for windows
t_loss_dur = duration(seconds(t_loss), 'Format','hh:mm');

%% Time axis
% Robust trimming (avoid undefined U)
t_ax_dur = duration(seconds(t_ax), 'Format','hh:mm');
valid_idx = find(P ~= 0, 1, 'first');
if isempty(valid_idx)
    valid_idx = 1;
end
t_lim = [t_ax_dur(valid_idx), t_ax_dur(end)];

%% Visualization

figure('Position',[200 100 1000 700])
t = tiledlayout(3,1,'TileSpacing','compact','Padding','compact');

% --- Power ---
ax(1) = nexttile;
% plot(t_ax, P, 'Color', [0 0.4470 0.7410 0.12]); hold on
plot(t_ax_dur, upperP, 'Color', [0 0.4470 0.7410], 'LineWidth', 2) 
ylabel('Power [W]')
title('Power')
grid on
% ylim([0.725 0.735])
% ylim(prctile(P,[1 99])) % robust scaling

% --- Temperature ---
ax(2) = nexttile;
plot(t_ax_dur, temperature, 'Color', [0.8500 0.3250 0.0980 0.05]); hold on
plot(t_ax_dur, tempSmooth, 'Color', [0.8500 0.3250 0.0980], 'LineWidth', 3)
ylabel('Temperature [°C]')
title('Implant Temperature')
grid on
ylim([27.5 30.5])


% --- Packet loss ---
ax(3) =nexttile;
% plot(t_ax, lostSamplePercentage, 'Color', [0 0.4470 0.7410 0.15]); hold on
% plot(t_ax, upperLossSmooth, 'Color', [0 0 0], 'LineWidth', 2)
plot(t_loss_dur,lostSamplePercentage,'Color', [0 0 0], 'LineWidth', 2)
yline(avgPL, 'LineWidth',1.2, 'LineStyle','--', 'Color','r')
ylabel('Packet-Loss [%]')
title(sprintf('Average Packet-Loss over %ds window', PL_twindow))
grid on

xlabel('Time [hours]')

% --- Link axes ---
linkaxes(ax,'x')
xlim(t_lim)

fontsize(12,'points')