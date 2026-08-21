% SUB_C04 - Subject-specific analysis script
%
% Author: Frederik Lampert
% Institution: Mayo Clinic
% Year: 2026
%
% This script configures and runs the BIC_data processing and analysis
% workflow for one subject. Update the placeholder paths below before use.
%
% Data availability: the source recordings are available in BCI2000 .dat
% format through OpenNeuro (doi:10.18112/openneuro.ds004624.v3.0.0).
%
% Requirements:
%   - MATLAB R2023b or later
%   - BIC_data class and associated project functions
%   - Any third-party utilities used below (for example progressbar and
%     custom colormap functions) must be added to the MATLAB path.
%
% NOTE: Subject-specific dates, channel assignments, excluded dates, and
% stimulation folders are retained from the original analysis and should
% be reviewed before adapting this script to another dataset.

clc; close all; clear all;
%% Configure paths
% Replace these placeholders with local paths before running the script.
% Keep project code, third-party dependencies, and data outside the script so
% the analysis remains portable acrosss systems.
projectPath = '/path/to/project/';
dependencyPath = '/path/to/dependencies/';
dataPath = '/path/to/data/';

% Project code containing BIC_data and associated helper functions.
addpath(genpath(projectPath));

% Third-party/helper functions used by the analysis (for example progressbar,
% legend utilities, spike detection code, and custom colormaps).
addpath(genpath(dependencyPath));

%% Select subject data directory
% Select the folder containing the subject's BCI2000 .dat recordings.
selpath = uigetdir(dataPath, 'Select subject data directory');
if isequal(selpath, 0)
    error('No data directory selected.');
end
%% Initialize the subject object
% If a previously saved subject object is present, it is reloaded. Otherwise,
% a new BIC_data object is initialized using the subject-specific metadata below.
if exist(fullfile(selpath,'Willie.mat'))
    temp = load(fullfile(selpath,'Willie.mat'));
    willie = temp.objCopy;
    willie = willie.loadData(selpath);
else
    willie = BIC_data(selpath);
    % Set subject-specific info
    willie.subjectName = 'Willie';
    willie.implantDate = datetime('Jul 28, 2025');
    willie.excludedDates = datetime({'Jul 28, 2025'});
    willie.knownBadChannels = [1:8];
    channelLeadMap = dictionary( ...
        'displaced', struct('indices', 1:8, 'labels',  "Ch" + (1:8)), ...
        'occipital',struct('indices', 9:16, 'labels',  "Ch" + (9:16)), ...
        'thalamus',struct('indices', 17:24, 'labels',  "Ch" + (17:24)),...
        'lobus-piriformis', struct('indices', 25:32, 'labels',  "Ch" + (25:32)) );
    willie.leadAssignment = channelLeadMap;
    willie.rerefSpecs.construct = 'combined';
    willie.rerefSpecs.specification = 'G-G-D-S';
    willie.rerefSpecs.ReRefCh = [9 17 25];  % example channel group
    willie.LongtermStimFolders = {'longterm-stim/LongtermStim50Hz_LANT001';...
        'longterm-stim/LongtermStim100Hz_LANT001';...
        'longterm-stim/LongtermStim125Hz_LANT001';...
        'longterm-stim/LongtermStim133Hz_LANT001';...
        'longterm-stim/LongtermStim200Hz_LANT001';...
        'longterm-stim/baseline001'};
    clear channelLeadMap
end
clear temp
%% Search for New data
% Search for all recordings
willie = willie.searchForNewData();
% Search for all Impedance files
willie = willie.searchForNewImp();
willie = willie.synchronizeFileStatus();
%% Process newly discovered recordings
% These steps may be computationally intensive and may overwrite/update cached
% processing results depending on the implementation of the BIC_data methods.
% Keep track of file processing status
willie = willie.synchronizeFileStatus();
tic; 
willie = willie.processRecordings(); % Initialize recordings structure and extract parameters from the files
% willie = willie.processLongtermStim();  % Initialize recordings structure and extract parameters from the files
willie = willie.synchronizeFileStatus();
willie = willie.extractSpectralFeatures([],1); toc;% Maybe optional input to specify the spectral calculation parameters?
willie = willie.getSpectralStats([],1);
toc;

%% Process long-term stim data
willie = willie.synchronizeFileStatus();
tic;
willie = willie.processLongtermStim();  
toc;

%% Detect spikes
% Spike detection depends on the configured detector and rereferencing methods.
% Review detector settings in BIC_data before reproducing this analysis.
willie = willie.synchronizeFileStatus();
willie = willie.detectSpikes([],true);

%% % Save the data
willie.saveData("object")

%% Summarize the recordings
[willie, sumTable] = willie.summarizeSubject();

%% Extract and all spectra and Impedance
allPSDs = willie.collectPSDs();
impTable = willie.summarizeImpedance();

%% Plot Impedance development
willie.plotImpedanceDevelopment(impTable)

%% Plot average PSD
% [avgPSD, PSDstd] =  willie.plotAveragePSD(allPSDs);
[avgPSD, PSDstd]  = willie.plotAveragePSD(allPSDs, true,'mean'); % Omit bad channels, use median instead of mean
ylim([-25 50])
axis square

%% Plot PSD evolution over time
% The custom plasma colormap is not part of base MATLAB. Replace it with a
% built-in colormap if the corresponding dependency is unavailable.
N = 50;
cmap = flipud(plasma(N+3)); %flipud(crameri('lapaz',N+2)); 
cmap(1:3,:) = []; %get rid of the white
willie.plotPSDoverTime(allPSDs, N, false, 'median', 'Native',cmap)
clear N

%% Plot RMS development
willie.plotStats('leadCAR','TotalRMS') % RMS of power spectra
% willie.plotStats('leadCAR','RMS')      % % RMS of time series
% xtickformat('MMM yyyy')
clim([0 500])

%% Plot Recording Capabilities
willie.plotRecordingCapabilities(impTable, allPSDs, 'leadCAR');

%% Plot packet loss count
f = figure('Name','Packet-Loss', 'Position', [100 100 500 800]);
bxplt = boxchart(getField(willie.recordings, 'PacketLoss'));
ylabel('Packet-loss [%]')
ylim([0 12])

%% Plots recording times
willie.plotRecordingTimes(willie.recordings)

%% Plot spike distribution
willie.plotSpikeDistribution('Native')

willie.plotNormalizedSpikeCounts('Native')
% seizure = datetime('01-Aug-2024 17:06:32', 'InputFormat','dd-MMM-uuuu HH:mm:ss');
% xline(seizure, '--', 'Color','red', 'LineWidth',2)
% legend({'Spike count','Data Trend', 'Seizure date'})

%% Plot spike examples
willie.plotSpikeExamples('leadCAR')

%% Plot Bad channel occurrence
badChDevelopment = willie.plotBadChannelOccurence();
xline(datetime({'Aug 15, 2025'; 'Jan 21, 2026';}), 'LineStyle','--', 'LineWidth',2, 'Color','r')
legend({'Number of functional channels','Max number of Channels', 'CT'})
filename = 'badChDevelopment.xlsx';
writetable(badChDevelopment, fullfile(selpath, filename))

%% #############################################
% Extract & Visualize long-term stimulation data
% ##############################################
%% Collect extracted features
[willie, stimFeatures] = willie.collectLongtermStimFeatures([]); %,'BP','Ch12-16'
%% Summarize in table the total recorded time for each type of stimulation
[stimSummary, stimTable] = willie.summarizeLongtermStim('BP');

%% Visualize long-term stim spectra
close all
refType = 'BP';
plotGroups = {'baseline'; ...
            'LongtermStim50Hz_LANT';...
            'LongtermStim100Hz_LANT';...
            'LongtermStim125Hz_LANT'};%;...
            % 'LongtermStim133Hz_LANT'};
plotData = willie.plotLongtermStim(refType,plotGroups);

%% Calculate statistics
mdl = fitglm(stimTable, ...
    'spikeCount ~ condition + recordingOrder + daysSinceImplant', ...
    'Distribution','poisson', ...
    'Offset',stimTable.logExposure);

%% % Plot additional boxplot for the spikeCount
    categs = categories(plotData.SpikeCount.cats);
    clrs = lines(numel(categs));
    
    figure('Position', [200 100 1000 500]);
    hold on
    
    for i = 1:numel(categs)
        thisMask = plotData.SpikeCount.cats == categs{i};
        boxchart(plotData.SpikeCount.cats(thisMask) , plotData.SpikeCount.vals(thisMask) , ...
            'BoxFaceColor', clrs(i,:), ...
            'MarkerStyle', 'none');
        
        xJitter = i + 0.075 * randn(sum(thisMask),1);
        scatter(xJitter, plotData.SpikeCount.vals(thisMask), ...
            35, ...
            clrs(i,:), ...
            'filled', ...
            'MarkerFaceAlpha', 0.3, ...
            'MarkerEdgeAlpha', 0.3);
    end
    
    ylabel('Spikes / (channel × min)')
    title('Normalized spike count')
    grid on
    
    ax = gca;
    ax.TickLabelInterpreter = 'none';
    fontsize(16, 'points');
%%  Plot additional boxplot for the spikeCount

% Keep only requested groups
stimTable2Plot = stimTable(ismember(string(stimTable.stimGroup), plotGroups), :);

% Remove unused categories and preserve desired order
stimTable2Plot.stimGroup = categorical( ...
    string(stimTable2Plot.stimGroup), ...
    plotGroups, ...
    'Ordinal', true);

categs = categories(stimTable2Plot.stimGroup);
clrs = lines(numel(categs));

figure('Position', [200 100 1000 500]);
hold on

for i = 1:numel(categs)
    thisMask = stimTable2Plot.stimGroup == categs{i};

    vals = stimTable2Plot.spikeCount(thisMask) ./ stimTable2Plot.logExposure(thisMask);
    cat = stimTable2Plot.stimGroup(thisMask);

    boxchart(cat, vals, ...
        'BoxFaceColor', clrs(i,:), ...
        'MarkerStyle', 'none');

    xJitter = i + 0.075*randn(sum(thisMask),1);
    scatter(xJitter, vals, ...
        100, ...
        clrs(i,:), ...
        'filled', ...
        'MarkerFaceAlpha', 0.3, ...
        'MarkerEdgeAlpha', 0.3);
end

ylabel('Spikes / log(channel × min)')
title('Normalized spike count')
grid on

ax = gca;
ax.TickLabelInterpreter = 'none';
fontsize(16,'points');
