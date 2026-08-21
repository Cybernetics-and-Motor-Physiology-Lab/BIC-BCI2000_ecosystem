% SUB_C00 - Subject-specific analysis script
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
if exist(fullfile(selpath,'Laika.mat'))
    temp = load(fullfile(selpath,'Laika.mat'));
    laika = temp.objCopy;
    laika = laika.loadData(selpath);
else
    laika = BIC_data(selpath);
    % Set subject-specific info
    laika.subjectName = 'Laika';
    laika.implantDate = datetime('April 23, 2021');
    laika.excludedDates = datetime({});
    channelLeadMap = dictionary( ...
        'left', struct('indices', 1:16, 'labels',  "Ch" + (1:16)), ...
        'right',struct('indices', 17:32, 'labels',  "Ch" + (17:32)));
    laika.leadAssignment = channelLeadMap;
    laika.rerefSpecs.construct = 'surface_wide';
    laika.rerefSpecs.specification = 'NN';
    laika.rerefSpecs.ReRefCh = [2 18];  
    laika.LongtermStimFolders = {};

end
clear temp
%% Search for New data
% Search for all recordings
laika = laika.searchForNewData();
% Search for all Impedance files
laika = laika.searchForNewImp();

%% Process newly discovered recordings
% These steps may be computationally intensive and may overwrite/update cached
% processing results depending on the implementation of the BIC_data methods.
tic; 
% Keep track of file processing status
laika = laika.synchronizeFileStatus();
laika = laika.processRecordings(); % Initialize recordings structure and extract parameters from the files

laika = laika.synchronizeFileStatus();
laika = laika.calculatePSDs(); % Maybe optional input to specify the spectral calculation parameters?
laika = laika.getSpectralStats();
toc;
laika = laika.synchronizeFileStatus();

%% Summarize the recordings
[laika,sumTable] = laika.summarizeSubject();

%% Save the data
laika.saveData()

%% Extract and all spectra and Impedance
allPSDs = laika.collectPSDs();

emptyFields = cellfun(@isempty, {allPSDs.name});
allPSDs(emptyFields) = [];

%% Plot average PSD
% [avgPSD, PSDstd] =  belka.plotAveragePSD(allPSDs);
[avgPSD, PSDstd]  = laika.plotAveragePSD(allPSDs, true,'mean'); % Omit bad channels, use median instead of mean
axis square

%% Plot PSD evolution over time
% The custom plasma colormap is not part of base MATLAB. Replace it with a
% built-in colormap if the corresponding dependency is unavailable.
% belka.plotPSDoverTime()
N = 20;
% cmap = cool(N);
cmap = flipud(plasma(N+3)); %flipud(crameri('lapaz',N+2)); 
cmap(1:3,:) = []; %get rid of the white
laika.plotPSDoverTime(allPSDs, N, true, 'median', 'leadCAR',cmap)
ylim([-25 10])
clear N

%% Plot packet loss count
f = figure('Name','Packet-Loss', 'Position', [100 100 1500 800]);
bxplt = boxchart(getField(laika.recordings, 'PacketLoss'));
ylabel('Packet-loss [%]')
ylim([0 10])
axis square

%% Plot RMS development
laika.plotStats('leadCAR','TotalRMS') %  plot RMS of power spectra
xticks(impTable.Date(1) : calmonths(2) : impTable.Date(end));
xtickformat('MMM yyyy')
clim([0 500])

%% Plots recording times
laika.plotRecordingTimes(laika.recordings)

%% Plot Bad channel occurrence
badChDevelopment = laika.plotBadChannelOccurence();
legend({'Number of functional channels','Max number of Channels', 'CT'})
filename = 'badChDevelopment.xlsx';
writetable(badChDevelopment, fullfile(selpath, filename))
