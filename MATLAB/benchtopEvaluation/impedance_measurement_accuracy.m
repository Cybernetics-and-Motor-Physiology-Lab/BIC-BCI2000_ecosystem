%% Impedance experiment – boxplot visualization
% MATLAB R2023b

clear; clc;

filename = "./BIC-BCI2000_ecosystem/data/impedanceMeasurmentAccuracy/Impedance_formated_for_Matlab.xlsx";

%% Load tables from each sheet
T_Rvals = readtable(filename, 'Sheet', 'Resistors');
T_Cvals = readtable(filename, 'Sheet', 'Capacitors');
T_R = readtable(filename, 'Sheet', 'R');
T_varR = readtable(filename, 'Sheet', 'varR');
T_varC = readtable(filename, 'Sheet', 'varC');

%% Labels from column names
RLabels = T_Rvals.Properties.VariableNames;
CLabels = T_Rvals.Properties.VariableNames;
varRLabels = T_varR.Properties.VariableNames;
varCLabels = T_varC.Properties.VariableNames;

%% Impedance meausrement parameters
Rpar = T_Rvals.ActualValue(5);
tr = 0.02e-3;      % 0.02 ms;

%% Create figure with tiled layout
figure('Position',[100 100 1200 500]);
t = tiledlayout(1,3, ...
    'TileSpacing','compact', ...
    'Padding','compact');

% 1) Pure resistor measurements
nexttile
boxchart(T_R.MeasuredValue); hold on;
plot([0.5 1.5], [T_R.ActualValue(1) T_R.ActualValue(1)], 'Color', [0.9294    0.6941    0.1255], 'LineStyle', '--', 'LineWidth', 1.5);
text(1.5, T_R.ActualValue(1) + 2, sprintf('Measured value: %.0f \\Omega',T_R.ActualValue(1)), ...
     'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom', ...
     'FontSize', 12);
title("Pure resistor measurements");
ylabel("Impedance [\Omega]");
grid on

% 2) RC impedance – varying resistance
C = T_Cvals.ActualValue(2); % For thi measurement C set  to 220nF
nexttile
boxchart(T_varR{:,2:5}); hold on;
for i=1:4
    Rseries = T_Rvals.ActualValue(i);
    [Zmag, Zphase_deg, feq] = calculateZ(C, Rseries, tr, Rpar);
    plot([i-0.5 i+0.5], [Zmag Zmag], 'Color', [0.9294    0.6941    0.1255], 'LineStyle', '--', 'LineWidth', 1.5);
    text(i+0.5, Zmag+2, sprintf('Calculated value: %.0f \\Omega',Zmag), ...
        'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom', ...
        'FontSize', 12);
end
title("RC impedance (varying resistance)");
ylabel("Impedance [\Omega]");
grid on

% 3) RC impedance – varying capacitance
nexttile
Rseries = T_Rvals.ActualValue(2); % For thi measurement R set  to 470Ω
boxchart(T_varC{:,2:4}); hold on;
for i=1:3
    C = T_Cvals.ActualValue(i);
    [Zmag, Zphase_deg, feq] = calculateZ(C, Rseries, tr, Rpar);
    plot([i-0.5 i+0.5], [Zmag Zmag], 'Color', [0.9294    0.6941    0.1255], 'LineStyle', '--', 'LineWidth', 1.5);
    text(i+0.5, Zmag+2, sprintf('Calculated value: %.0f \\Omega',Zmag), ...
        'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom', ...
        'FontSize', 12);
end
title("RC impedance (varying capacitance)");
ylabel("Impedance [\Omega]");
grid on

% Global title
title(t, "Impedance measurements across experimental conditions");

