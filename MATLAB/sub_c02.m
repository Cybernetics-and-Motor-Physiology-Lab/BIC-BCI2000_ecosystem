% SUB_C02 - Subject-specific analysis script
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
if exist(fullfile(selpath,'Billy.mat'))
    temp = load(fullfile(selpath,'Billy.mat'));
    billy = temp.objCopy;
    billy = billy.loadData(selpath);
else
    billy = BIC_data(selpath);
    % Set subject-specific info
    billy.subjectName = 'Billy';
    billy.implantDate = datetime('May 13, 2024');
    billy.excludedDates = datetime({'May 13, 2024','Jun 03, 2024','Aug 01, 2024'});
    billy.knownBadChannels = [1 2 5 6];
    channelLeadMap = dictionary( ...
        'LANT', struct('indices', 1:8, 'labels',  "Ch" + (1:8)), ...
        'LHPC',struct('indices', 9:16, 'labels',  "Ch" + (9:16)), ...
        'RANT',struct('indices', 17:24, 'labels',  "Ch" + (17:24)),...
        'RHPC', struct('indices', 25:32, 'labels',  "Ch" + (25:32)) );
    billy.leadAssignment = channelLeadMap;
    billy.rerefSpecs.construct = 'depth';
    billy.rerefSpecs.specification = 'D-S-D-S';
    billy.rerefSpecs.ReRefCh = [8 16 24 25];  % example channel group
    billy.LongtermStimFolders = {'longterm-stim/LongtermStim001'};
end
clear temp channelLeadMap
%% Search for New data
% Search for all recordings
billy = billy.searchForNewData();
% Search for all Impedance files
billy = billy.searchForNewImp();
billy = billy.synchronizeFileStatus();
%% Process newly discovered recordings
% These steps may be computationally intensive and may overwrite/update cached
% processing results depending on the implementation of the BIC_data methods.
% Keep track of file processing status
billy = billy.synchronizeFileStatus();
tic; 
billy = billy.processRecordings(); % Initialize recordings structure and extract parameters from the files
billy = billy.processLongtermStim();  % Initialize recordings structure and extract parameters from the files
billy = billy.synchronizeFileStatus();
billy = billy.extractSpectralFeatures();% Maybe optional input to specify the spectral calculation parameters?
billy = billy.getSpectralStats();
toc;

% Process long-term stim
tic
billy = billy.synchronizeFileStatus();
billy = billy.processLongtermStim();  % Initialize recordings structure and extract parameters from the files
toc;

% Detect spikes
billy = billy.synchronizeFileStatus();
billy = billy.detectSpikes([],true);

% Save the data
billy.saveData()

% Summarize the recordings
[billy, sumTable] = billy.summarizeSubject();

%% Extract and all spectra and Impedance
allPSDs = billy.collectPSDs();
impTable = billy.summarizeImpedance();

%% Plot Impedance development
billy.plotImpedanceDevelopment(impTable)

%% Plot average PSD
% [avgPSD, PSDstd] =  billy.plotAveragePSD(allPSDs);
[avgPSD, PSDstd]  = billy.plotAveragePSD(allPSDs, true,'mean'); % Omit bad channels, use median instead of mean
ylim([-25 50])

%% Plot PSD evolution over time
% The custom plasma colormap is not part of base MATLAB. Replace it with a
% built-in colormap if the corresponding dependency is unavailable.
N = 30;
cmap = flipud(plasma(N+3)); %flipud(crameri('lapaz',N+2)); 
cmap(1:3,:) = []; %get rid of the white
billy.plotPSDoverTime(allPSDs, N, true, 'median', 'Native',cmap)
ylim([-25 10])
clear N

%% Plot packet loss count
f = figure('Name','Packet-Loss', 'Position', [100 100 500 800]);
bxplt = boxchart(getField(billy.recordings, 'PacketLoss'));
ylabel('Packet-loss [%]')
ylim([0 10])

%% Plot RMS development
billy.plotStats('leadCAR','TotalRMS') % RMS of power spectra
% billy.plotStats('leadCAR','RMS')      % % RMS of time series
xtickformat('MMM yyyy')
clim([0 500])

%% Plot Recording Capabilities
billy.plotRecordingCapabilities(impTable, allPSDs, 'leadCAR');

%% Plot avg PSD over time
% billy.plotPSDoverTime()
billy.plotPSDoverTime('leadCAR','median'); clim([-20 20]); ylim([0 250])

%% Extract and plot long-term stimulation spectra
% This section requires long-term stimulation files and any intermediate
% variables referenced below to have been generated in earlier processing steps.
longtermPSDs = billy.collectLongtermStimPSDs();

refType = 'leadCAR';
billy.plotLongtermStim(avgPSD.(refType), longtermPSDs, true, refType)

%% Plots recording times
billy.plotRecordingTimes(billy.recordings)

%% Plot spike distribution
billy.plotSpikeDistribution('BP')

billy.plotNormalizedSpikeCounts('BP')
seizure = datetime('01-Aug-2024 17:06:32', 'InputFormat','dd-MMM-uuuu HH:mm:ss');
xline(seizure, '--', 'Color','red', 'LineWidth',2)
legend({'Spike count','Data Trend', 'Seizure date'})


%% Plot spike examples
billy.plotSpikeExamples('BP')

%% Plot Bad channel occurrence
badChDevelopment = billy.plotBadChannelOccurence();
xline(datetime({'May 13, 2024'}), 'LineStyle','--', 'LineWidth',2, 'Color','r')
legend({'Number of functional channels','Max number of Channels', 'CT'})
filename = 'badChDevelopment.xlsx';
writetable(badChDevelopment, fullfile(selpath, filename))
