% SUB_S05 - Subject-specific analysis script
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
if exist(fullfile(selpath,'Mushka.mat'))
    temp = load(fullfile(selpath,'Mushka.mat'));
    mushka = temp.objCopy;
    mushka = mushka.loadData(selpath);
else
    mushka = BIC_data(selpath);
    % Set subject-specific info
    mushka.subjectName = 'Mushka';
    mushka.implantDate = datetime('January 12, 2026');
    mushka.excludedDates = datetime({'January 26, 2026', 'February 13,2026'});
    channelLeadMap = dictionary( ...
        'frontal', struct('indices', 1:8, 'labels',  "Ch" + (1:8)), ...
        'temporal',struct('indices', 9:16, 'labels',  "Ch" + (9:16)), ...
        'occipital',struct('indices', 17:24, 'labels',  "Ch" + (17:24)),...
        'thalamus', struct('indices', 25:32, 'labels',  "Ch" + (25:32)) );
    mushka.leadAssignment = channelLeadMap;
    mushka.rerefSpecs.construct = 'combined';
    mushka.rerefSpecs.specification = 'G-G-G-D';
    mushka.rerefSpecs.ReRefCh = [1 9 17 25];  % example channel group
    mushka.LongtermStimFolders = {};

end
clear temp
%% Search for New data
% Search for all recordings
mushka = mushka.searchForNewData();
% Search for all Impedance files
mushka = mushka.searchForNewImp();

%% Process newly discovered recordings
% These steps may be computationally intensive and may overwrite/update cached
% processing results depending on the implementation of the BIC_data methods.
tic; 
% Keep track of file processing status
mushka = mushka.synchronizeFileStatus();
mushka = mushka.processRecordings(); % Initialize recordings structure and extract parameters from the files

mushka = mushka.synchronizeFileStatus();
mushka = mushka.extractSpectralFeatures(); % Maybe optional input to specify the spectral calculation parameters?
mushka = mushka.getSpectralStats([],1);
toc;
mushka = mushka.synchronizeFileStatus();

%% Detect spikes
% Spike detection depends on the configured detector and rereferencing methods.
% Review detector settings in BIC_data before reproducing this analysis.
mushka = mushka.synchronizeFileStatus();
% belka = belka.detectSpikes('BP',true);
mushka = mushka.detectSpikes([],true);

%% Summarize the recordings
[mushka,sumTable] = mushka.summarizeSubject();

%% Save the data
mushka.saveData()

%% Extract and all spectra and Impedance
allPSDs = mushka.collectPSDs();
impTable = mushka.summarizeImpedance();

emptyFields = cellfun(@isempty, {allPSDs.name});
allPSDs(emptyFields) = [];
%% Plot Impedance development
mushka.plotImpedanceDevelopment(impTable)

%% Plot average PSD
% [avgPSD, PSDstd] =  belka.plotAveragePSD(allPSDs);
[avgPSD, PSDstd]  = mushka.plotAveragePSD(allPSDs, true,'mean'); % Omit bad channels, use median instead of mean
axis square

%% Plot PSD evolution over time
% The custom plasma colormap is not part of base MATLAB. Replace it with a
% built-in colormap if the corresponding dependency is unavailable.
% belka.plotPSDoverTime()
N = 50;
% cmap = cool(N);
cmap = flipud(plasma(N+1)); %flipud(crameri('lapaz',N+2)); 
cmap(1,:) = []; %get rid of the most bright color
mushka.plotPSDoverTime(allPSDs, N, false, 'median', 'Native',cmap)
% belka.plotPSDoverTime('leadCAR','median'); clim([-20 20]);

clear N

%% Plot packet loss count
f = figure('Name','Packet-Loss', 'Position', [100 100 1500 800]);
bxplt = boxchart(getField(mushka.recordings, 'PacketLoss'));
ylabel('Packet-loss [%]')
ylim([0 10])
axis square

%% Plot RMS development
mushka.plotStats('leadCAR','TotalRMS') %  plot RMS of power spectra
% belka.plotStats('leadCAR','RMS')      % % RMS of time series
% xticks(impTable.Date(1) : calmonths(2) : impTable.Date(end));
xtickformat('dd/MMM yyyy');
clim([0 500])
%% Plot Recording Capabilities
mushka.plotRecordingCapabilities(impTable, allPSDs, 'leadCAR');

%% Extract and plot long-term stimulation spectra
% This section requires long-term stimulation files and any intermediate
% variables referenced below to have been generated in earlier processing steps.
% longtermPSDs = belka.collectLongtermStimPSDs();

refType = 'leadCAR';
mushka.plotLongtermStim(avgPSD.(refType), longtermPSDs, true, refType)

%% Plots recording times
mushka.plotRecordingTimes(mushka.recordings)

%% Plot spike distribution
mushka.plotSpikeDistribution('BP')

mushka.plotNormalizedSpikeCounts('BP')
xticks(impTable.Date(1) : calmonths(4) : impTable.Date(end));
xtickformat('MMM yyyy')
%% Plot Bad channel occurrence
badChDevelopment = mushka.plotBadChannelOccurence();
xline(datetime({'Jan 12, 2026'; 'Apr 06, 2026'}), 'LineStyle','--', 'LineWidth',2, 'Color','r')
legend({'Number of functional channels','Max number of Channels', 'CT'})
filename = 'badChDevelopment.xlsx';
writetable(badChDevelopment, fullfile(selpath, filename))
