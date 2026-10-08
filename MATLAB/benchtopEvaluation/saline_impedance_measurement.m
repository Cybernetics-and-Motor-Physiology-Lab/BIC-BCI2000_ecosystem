%% Impedance experiment – boxplot visualization
% MATLAB R2023b
clear; clc;

%% Load data
fileName = "../../data/saline-ImpedanceMeasurments/Impedance_measurments_saline.xlsx";
raw = readcell(fileName, 'Sheet', 'Sheet1');

% Row 1:
% Column B:O = NaCl concentrations
% Column P   = distilled water (0% NaCl)
%
% Row 2 contains condition names
% Rows 3:end contain measurements

concentration = cell2mat(raw(1, 2:15));     % C1-C14
impedance     = cell2mat(raw(3:end-2, 2:15)); % saline measurements

%% Distilled water measurements
waterRaw = raw(3:end-2, 16);

% Convert "Inf" stored as text to numerical Inf
waterImpedance = nan(size(waterRaw));

for i = 1:numel(waterRaw)
    if isnumeric(waterRaw{i})
        waterImpedance(i) = waterRaw{i};
    elseif ischar(waterRaw{i}) || isstring(waterRaw{i})
        if strcmpi(string(waterRaw{i}), "Inf")
            waterImpedance(i) = Inf;
        else
            waterImpedance(i) = str2double(string(waterRaw{i}));
        end
    end

end

%% Replace Inf with 20 kOhm
impedance(isinf(impedance)) = 20000;
waterImpedance(isinf(waterImpedance)) = 20000;

%% Calculate statistics
meanImpedance = mean(impedance, 1, 'omitnan');
stdImpedance  = std(impedance, 0, 1, 'omitnan');

%% Sort concentrations from low -> high
% Makes the plot easier to interpret from left to right
[concentration, sortIdx] = sort(concentration);

impedance     = impedance(:, sortIdx);
meanImpedance = meanImpedance(sortIdx);
stdImpedance  = stdImpedance(sortIdx);

%% Position for distilled water
% 0 cannot be displayed on a logarithmic axis.
% Place distilled water one logarithmic step to the left of the
% lowest measured NaCl concentration.

waterX = min(concentration) / 2.5;

%% Create figure
fig = figure( ...
    'Color', 'w', ...
    'Position', [100 100 1100 700]);

ax = axes(fig);
hold(ax, 'on');

%% Mean +/- standard deviation
errorbar(ax, ...
    concentration, ...
    meanImpedance, ...
    stdImpedance, ...
    '_', ...
    'LineStyle', 'none', ...
    'MarkerFaceColor', 'none', ...
    'MarkerEdgeColor', [0.20 0.45 0.70], ...
    'Color', [0.20 0.45 0.70], ...
    'MarkerSize', 12, ...
    'LineWidth', 1.2, ...
    'CapSize', 8);

%% Distilled water
scatter(ax, ...
    waterX, ...
    20000, ...
    65, ...
    '^', ...
    'k', ...
    'filled');

text(ax, waterX*1.05, 20000, ...
    '  >20 k\Omega', ...
    'FontSize', 11, ...
    'VerticalAlignment', 'middle');

%% Axes
ax.XScale = 'log';
ax.YScale = 'log';

% Explicit concentration ticks
ax.XTick = [waterX concentration];
ax.XTickLabel = ["DI water", compose('%.3f', concentration)];
ax.XTickLabelRotation = 45;

xlabel(ax, 'NaCl concentration (%)');
ylabel(ax, 'Impedance (\Omega)');

%% Clean grid
grid(ax, 'on');
ax.XMinorGrid = 'off';
ax.YMinorGrid = 'off';

% Make major grid subtle
ax.GridAlpha = 0.15;
ax.GridLineStyle = '-';

% %% Tick marks
% ax.TickDir = 'out';
% ax.TickLength = [0.008 0.008];

%% General appearance
ax.FontSize = 11;
ax.LineWidth = 1;
ax.Box = 'off';

% Remove title for manuscript figure
title(ax, '');

%% Limits
xlim(ax, [waterX/1.3 max(concentration)*1.25]);
ylim(ax, [100 25000]);

%% Figure
fig.Color = 'w';
%% Optional title
title(ax, 'BIC impedance as a function of NaCl concentration');