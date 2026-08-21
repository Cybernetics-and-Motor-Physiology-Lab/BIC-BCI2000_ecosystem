% SUB_C03 - Subject-specific analysis script
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
if exist(fullfile(selpath,'Strelka.mat'),'file')
    temp = load(fullfile(selpath,'Strelka.mat'));
    strelka = temp.objCopy;
    strelka = strelka.loadData(selpath);
else
    strelka = BIC_data(selpath);
    % Set subject-specific info
    strelka.subjectName = 'Strelka';
    strelka.implantDate = datetime('Nov 04, 2024');
    strelka.excludedDates = datetime({'Nov 04, 2024','Nov 22, 2024'});
    channelLeadMap = dictionary( ...
        'frontal', struct('indices', 1:12, 'labels',  "Ch" + (1:12)), ...
        'temporal',struct('indices', 13:22, 'labels',  "Ch" + (13:22)), ...
        'occipital',struct('indices', 23:32, 'labels',  "Ch" + (23:32)));
    strelka.leadAssignment = channelLeadMap;
    strelka.rerefSpecs.construct = 'surface';
    strelka.rerefSpecs.specification = 'NN';
    strelka.rerefSpecs.ReRefCh = [4 22 28];  % example channel group
    strelka.LongtermStimFolders = {'longterm/Stim_Grid1001';...
        'longterm/Stim_Grid2001';...
        'longterm/Stim_Grid3001'};
end
clear temp
%% Search for New data
% Search for all recordings
strelka = strelka.searchForNewData();
% Search for all Impedance files
strelka = strelka.searchForNewImp();

%% Process newly discovered recordings
% These steps may be computationally intensive and may overwrite/update cached
% processing results depending on the implementation of the BIC_data methods.
% Keep track of file processing status
strelka = strelka.synchronizeFileStatus();
tic; 
strelka = strelka.processRecordings(); % Initialize recordings structure and extract parameters from the files
strelka = strelka.processLongtermStim();  % Initialize recordings structure and extract parameters from the files
strelka = strelka.synchronizeFileStatus();
% strelka = strelka.calculatePSDs(); % Maybe optional input to specify the spectral calculation parameters?
strelka = strelka.getSpectralStats();
toc;

%% Detect spikes
% Spike detection depends on the configured detector and rereferencing methods.
% Review detector settings in BIC_data before reproducing this analysis.
strelka = strelka.synchronizeFileStatus();
% strelka = strelka.detectSpikes('BP',true);
strelka = strelka.detectSpikes([],true);
%% Summarize the recordings
[strelka,sumTable] = strelka.summarizeSubject();

%% Save the data
strelka.saveData()

%% Extract and all spectra and Impedance
allPSDs = strelka.collectPSDs();
impTable = strelka.summarizeImpedance();

%% Plot Impedance development
strelka.plotImpedanceDevelopment(impTable)

%% Plot average PSD
% [avgPSD, PSDstd] =  strelka.plotAveragePSD(allPSDs);
[avgPSD, PSDstd]  = strelka.plotAveragePSD(allPSDs, true,'mean'); % Omit bad channels, use median instead of mean
ylim([-25 50])
axis square

%% Plot PSD evolution over time
% The custom plasma colormap is not part of base MATLAB. Replace it with a
% built-in colormap if the corresponding dependency is unavailable.
N = 50;
cmap = flipud(plasma(N+3)); %flipud(crameri('lapaz',N+2)); 
cmap(1:3,:) = []; %get rid of the white
strelka.plotPSDoverTime(allPSDs, N, false, 'median', 'Native',cmap)
% ylim([-25 10])
clear N

%% Plot packet loss count
f = figure('Name','Packet-Loss', 'Position', [100 100 500 800]);
bxplt = boxchart(getField(strelka.recordings, 'PacketLoss'));
ylabel('Packet-loss [%]')
ylim([0 10])

%% Plot Recording Capabilities
strelka.plotRecordingCapabilities(impTable, allPSDs, 'leadCAR');

%% Plot RMS development
strelka.plotStats('leadCAR','TotalRMS') %  plot RMS of power spectra
% strelka.plotStats('leadCAR','STD')      % % RMS of time series
clim([0 500])

%% Extract and plot long-term stimulation spectra
% This section requires long-term stimulation files and any intermediate
% variables referenced below to have been generated in earlier processing steps.
% longtermPSDs = strelka.collectLongtermStimPSDs();

refType = 'leadCAR';
strelka.plotLongtermStim(avgPSD.(refType), longtermPSDs, true, refType)

%% Plots recording times
strelka.plotRecordingTimes(strelka.recordings)

%% Plot spike distribution
strelka.plotSpikeDistribution('BP')

strelka.plotNormalizedSpikeCounts('BP')

%% Plot Bad channel occurrence
badChDevelopment = strelka.plotBadChannelOccurence();
xline(datetime({'Nov 4, 2024'; 'Mar 03, 2025'; 'Jan 21, 2026';}), 'LineStyle','--', 'LineWidth',2, 'Color','r')
legend({'Number of functional channels','Max number of Channels', 'CT'})
filename = 'badChDevelopment.xlsx';
writetable(badChDevelopment, fullfile(selpath, filename))
