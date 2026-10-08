% LONGTERM_INVIVO_REC_CAPABILITY - longterm inVivo rec capability
%
% Author: Frederik Lampert
% Institution: Mayo Clinic
% Year: 2026
%
% This script reproduces an in vivo analysis workflow used for evaluation
% of long-term neural recording performance and signal evolution.
%
% Data availability: source recordings are available in BCI2000 .dat format
% through Dandi Archive (https://dandiarchive.org/dandiset/000571).
%
% Before running:
%   1. Replace placeholder paths with local project/data locations.
%   2. Add the BIC_data class and required helper functions to the MATLAB path.
%   3. Review subject-specific selections, exclusions, and plotting settings.
%
% MATLAB version: R2023b

clear all; close all; clc

%% Configure project path
% Replace with the local repository or source-code folder.
projectPath = '/path/to/project/';
addpath(genpath(projectPath));

%% Specify subjects and data path
dataPath = 'path/to/BIC-BCI2000_ecosystem/data';
% Subject/file selections below are analysis-specific; review before reuse.

subjects = {
    'sub-c01';
    'sub-c02';
    'sub-c03';
    'sub-c04';
    'sub-c05';
    };

subColors = lines(numel(subjects));

%% Load all packet loss files
PL = struct();

for i = 1:numel(subjects)
    sub = subjects{i};
    plPath = fullfile(dataPath, sub, sprintf('%s_PacketLoss.mat', sub));
    tmp = load(plPath);
    vnames = fieldnames(tmp);
    % Load first variable from the .mat file
    PL.(matlab.lang.makeValidName(sub)) = tmp.(vnames{1});
end

%% Plot results
% Plotting parameters below reproduce the analysis figures and may be adapted as needed. packet loss boxcharts with overlaid scatters
plVals = [];
plSubj = strings(0,1);

for i = 1:numel(subjects)
    sub = subjects{i};
    subField = matlab.lang.makeValidName(sub);
    vals = PL.(subField);
    vals = vals(:);
    plVals = [plVals; vals];
    plSubj = [plSubj; repmat(string(sub), numel(vals), 1)];
end

plCats = categorical(plSubj, string(subjects), 'Ordinal', true);

figure('Position', [200 100 1000 500]);
hold on

for i = 1:numel(subjects)
    thisMask = plCats == subjects{i};
    boxchart(plCats(thisMask), plVals(thisMask), ...
        'BoxFaceColor', subColors(i,:), ...
        'MarkerStyle', 'none');
    xJitter = i + 0.075 * randn(sum(thisMask),1);
    scatter(xJitter, plVals(thisMask), ...
        35, ...
        subColors(i,:), ...
        'filled', ...
        'MarkerFaceAlpha', 0.3, ...
        'MarkerEdgeAlpha', 0.3);
end

ylim([0 15])
ylabel('Packet loss [%]')
title('Packet loss across subjects')
grid on

ax = gca;
ax.TickLabelInterpreter = 'none';
fontsize(16, 'points');

%% Load bad channel tables and aggregate
plotData = struct();
plotData.Subject = strings(0,1);
plotData.Month = [];
plotData.TotalGoodCh = [];

for i = 1:numel(subjects)
    sub = subjects{i};
    tblPath = fullfile(dataPath, sub, 'badChDevelopment.xlsx');
    badChnlsT = readtable(tblPath);

    subField = matlab.lang.makeValidName(sub);
    plotData.(subField).Table = badChnlsT;
    plotData.(subField).Month = badChnlsT.Month;
    plotData.(subField).TotalGoodCh = badChnlsT.TotalGoodCh;
    plotData.Subject = [plotData.Subject; repmat(string(sub), height(badChnlsT), 1)];
    plotData.Month = [plotData.Month; badChnlsT.Month];
    plotData.TotalGoodCh = [plotData.TotalGoodCh; badChnlsT.TotalGoodCh];
end

%% Plot total good channels over time
figure('Position', [200 100 1000 500]);
hold on

for i = 1:numel(subjects)
    sub = subjects{i};
    subMask = plotData.Subject == string(sub);
    stairs(plotData.Month(subMask), plotData.TotalGoodCh(subMask), ...
        'Color', subColors(i,:), ...
        'LineWidth', 1.5, ...
        'DisplayName', sub);
end

% stairs(plotData.Month,plotData.TotalGoodCh,'LineWidth',2,'Marker','*')
yline(32,'LineStyle','--')
ylim([15 33])
            
xlabel('Month')
xlim('tight')
% xticks(plotData.Month(1) : calmonths(2) : plotData.Month(end));
xtickangle(45)
ylabel('Total good channels')
title('Good-channel development across subjects')
legend('Location', 'best', 'Interpreter', 'none')
grid on
fontsize(14, 'points');
