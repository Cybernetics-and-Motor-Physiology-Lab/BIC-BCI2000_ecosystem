% SUB_C01 - Subject-specific analysis script
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
if exist(fullfile(selpath,'Belka.mat'))
    temp = load(fullfile(selpath,'Belka.mat'));
    belka = temp.objCopy;
    belka = belka.loadData(selpath);
else
    belka = BIC_data(selpath);
    % Set subject-specific info
    belka.subjectName = 'Belka';
    belka.implantDate = datetime('March 20, 2023');
    belka.excludedDates = datetime({'Jun 02, 2023'; 'Jun 06, 2023'; 'Oct 31, 2023'; 'Nov 6, 2023'; 'Nov 8, 2023'; 'Jan 2, 2024';'Jan 15, 2024'; 'Jan 17, 2024'; 'Jan 18, 2024'; 'Jan 23, 2024'; 'Feb 26, 2024'; 'Feb 29, 2024';'May 02, 2025'});
    channelLeadMap = dictionary( ...
        'frontal', struct('indices', 1:12, 'labels',  "Ch" + (1:12)), ...
        'temporal',struct('indices', 13:22, 'labels',  "Ch" + (13:22)), ...
        'occipital',struct('indices', 23:32, 'labels',  "Ch" + (23:32)));
    belka.leadAssignment = channelLeadMap;
    belka.rerefSpecs.construct = 'surface';
    belka.rerefSpecs.specification = 'NN';
    belka.rerefSpecs.ReRefCh = [2 18 31];  
    belka.LongtermStimFolders = {'longterm/stim_recordings/Stim_6-12001';...
        'longterm/stim_recordings/Stim_13-18001';...
        'longterm/stim_recordings/Stim_23-28001'};

end
clear temp
%% Search for New data
% Search for all recordings
belka = belka.searchForNewData();
% Search for all Impedance files
belka = belka.searchForNewImp();

%% Process newly discovered recordings
% These steps may be computationally intensive and may overwrite/update cached
% processing results depending on the implementation of the BIC_data methods.
tic; 
% Keep track of file processing status
belka = belka.synchronizeFileStatus();
belka = belka.processRecordings(); % Initialize recordings structure and extract parameters from the files

belka = belka.synchronizeFileStatus();
belka = belka.extractSpectralFeatures(); % Maybe optional input to specify the spectral calculation parameters?
belka = belka.getSpectralStats();
toc;
belka = belka.synchronizeFileStatus();

%% Process long-term stim data
belka = belka.processLongtermStim();  

%% Detect spikes
% Spike detection depends on the configured detector and rereferencing methods.
% Review detector settings in BIC_data before reproducing this analysis.
belka = belka.synchronizeFileStatus();
% belka = belka.detectSpikes('BP',true);
belka = belka.detectSpikes([],true);

%% Summarize the recordings
[belka,sumTable] = belka.summarizeSubject();

%% Save the data
belka.saveData()

%% Extract and all spectra and Impedance
allPSDs = belka.collectPSDs();
impTable = belka.summarizeImpedance();

emptyFields = cellfun(@isempty, {allPSDs.name});
allPSDs(emptyFields) = [];
%% Plot Impedance development
belka.plotImpedanceDevelopment(impTable)
% set(gca, 'YDir', 'normal');

%% Plot average PSD
% [avgPSD, PSDstd] =  belka.plotAveragePSD(allPSDs);
[avgPSD, PSDstd]  = belka.plotAveragePSD(allPSDs, true,'mean'); % Omit bad channels, use median instead of mean
axis square

%% Plot PSD evolution over time
% The custom plasma colormap is not part of base MATLAB. Replace it with a
% built-in colormap if the corresponding dependency is unavailable.
% belka.plotPSDoverTime()
N = 50;
% cmap = cool(N);
cmap = flipud(plasma(N+3)); %flipud(crameri('lapaz',N+2)); 
cmap(1:3,:) = []; %get rid of the white
belka.plotPSDoverTime(allPSDs, N, false, 'median', 'Native',cmap)
% belka.plotPSDoverTime('leadCAR','median'); clim([-20 20]); 
% ylim([-25 10])
clear N

%% Plot packet loss count
f = figure('Name','Packet-Loss', 'Position', [100 100 1500 800]);
bxplt = boxchart(getField(belka.recordings, 'PacketLoss'));
ylabel('Packet-loss [%]')
ylim([0 10])
axis square

%% Plot RMS development
belka.plotStats('leadCAR','TotalRMS') %  plot RMS of power spectra
% belka.plotStats('Native','ThetaRMS')      % % RMS of time series
xticks(impTable.Date(1) : calmonths(2) : impTable.Date(end));
xtickformat('MMM yyyy')
clim([0 500])
%% Plot Recording Capabilities
belka.plotRecordingCapabilities(impTable, allPSDs, 'leadCAR');

%% Extract and plot long-term stimulation spectra
% This section requires long-term stimulation files and any intermediate
% variables referenced below to have been generated in earlier processing steps.
% longtermPSDs = belka.collectLongtermStimPSDs();

refType = 'leadCAR';
belka.plotLongtermStim(avgPSD.(refType), longtermPSDs, true, refType)

%% Plots recording times
belka.plotRecordingTimes(belka.recordings)

%% Plot spike distribution
belka.plotSpikeDistribution('BP')

belka.plotNormalizedSpikeCounts('BP')
xticks(impTable.Date(1) : calmonths(4) : impTable.Date(end));
xtickformat('MMM yyyy')
%% Plot Bad channel occurrence
badChDevelopment = belka.plotBadChannelOccurence();
xline(datetime({'Oct 24, 2024'; 'Mar 03, 2025'; 'Jan 11, 2026';}), 'LineStyle','--', 'LineWidth',2, 'Color','r')
legend({'Number of functional channels','Max number of Channels', 'CT'})
filename = 'badChDevelopment.xlsx';
writetable(badChDevelopment, fullfile(selpath, filename))
