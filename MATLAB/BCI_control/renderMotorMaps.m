%% renderMotorMaps
% Render cortical and MRI-slice motor maps for subject IMG.
%
% Authors: Jordan Bilderbeek, Frederik Lampert
% Institution: Mayo Clinic
% Year: 2026
%
% DESCRIPTION
%   This script visualizes motor-mapping results in two complementary views:
%     1. A 3-D cortical surface with electrode markers scaled and colored by
%        signed r^2 values.
%     2. Sagittal, coronal, and axial MRI slices with weighted electrode
%        overlays using mnl_seegview.
%
% DATA AVAILABILITY
%   Source recordings are available in BCI2000 .dat format through DandiArchive:
%   https://dandiarchive.org/dandiset/000571
%
% REQUIREMENTS
%   - MATLAB R2023b
%   - mnl_seegview:
%       https://github.com/MultimodalNeuroimagingLab/mnl_seegview
%   - mnl_ieegBasics:
%       https://github.com/MultimodalNeuroimagingLab/mnl_ieegBasics
%   - slanCM colormaps:
%       https://www.mathworks.com/matlabcentral/fileexchange/120088-200-colormaps
%   - Project-specific helper functions used below, including locs_DPRR,
%     ieeg_RenderGifti, ieeg_elAdd, ieeg_viewLight, seegview_sliceplot,
%     sv_weight_add, and kjm_printfig.
%
% NOTES
%   Electrode selections and task-specific result fields are subject-specific.
%   Review the configuration section before adapting this script to another
%   participant or motor-mapping condition.

clear; close all; clc

%% Configure paths
addpath('./functions/');
addpath(genpath('.../mnl_seegview')); % https://github.com/MultimodalNeuroimagingLab/mnl_seegview
addpath(genpath('path/to/mnl_ieegBasics/functions')); % https://github.com/MultimodalNeuroimagingLab/mnl_ieegBasics
addpath(genpath('.../slanCM')); % https://www.mathworks.com/matlabcentral/fileexchange/120088-200-colormaps 

%% Configure subject and analysis
dataPath  = '../data';
subjectID = 'IMG';
taskName = 'overtBCI_d3s2';

% Electrode metadata, motor-mapping results, and anatomical data.
electrodesTsv = fullfile(dataPath, 'BCI', 'IMG_electrodes.tsv');
rvalsFile     = fullfile(dataPath, 'BCI', 'IMG_r_mot_rvals_noBon.mat');
brainFile     = fullfile(dataPath, 'brainMAT', 'sub-IMG_brain.mat');

% Subject-specific subset of electrodes associated with the cortical array.
% Indices refer to the subset of electrodes that were used later for BIC recording.
% d3s2 (overt BCI) 74:101 | d4s2 (imagery BCI) 27:54
cortecDPIdx = 74:101; 

% Motor-mapping variable to visualize on the cortical surface.
r2Field = 'r_foot_HFB';

% 3-D electrode-marker size range.
minElSize = 5;
maxElSize = 20;

% MRI-slice visualization settings.
side = 'B';
sliceThickness = 5;      % mm
mrClims = [0 1];         % MRI display limits
pThresh = 0.05;          % significance threshold
threshFrac = 0.01;       % suppress weights below this fraction of maximum
saveSliceFigures = false;

%% Load electrode locations and analysis data
assert(isfile(electrodesTsv), 'Electrode file not found: %s', electrodesTsv);
assert(isfile(rvalsFile), 'Motor-mapping results not found: %s', rvalsFile);
assert(isfile(brainFile), 'Brain anatomy file not found: %s', brainFile);

load(rvalsFile);
load(brainFile);

%% Prepare electrode locations
locInfo = ieeg_readtableRmHyphens(electrodesTsv);
locs = [locInfo.x, locInfo.y, locInfo.z];

% Project electrode locations for slice-based visualization.
[dpChannels, dpLocs] = locs_DPRR(locs, 5, 'y');
locInfoDP = locInfo(dpChannels, :);

% Split the projected electrode table into the cortical-array subset and all
% remaining electrodes.
locInfoCortec = locInfoDP(cortecDPIdx, :);
locInfoRest = locInfoDP;
locInfoRest(cortecDPIdx, :) = [];


%% Render 3-D cortical motor map
fig = figure( ...
    'Position', [100 100 2000 1500], ...
    'Visible', 'on', ...
    'PaperUnits', 'inches', ...
    'PaperPosition', [0 0 20 15], ...
    'PaperPositionMode', 'auto', ...
    'PaperSize', [20 15], ...
    'PaperType', '<custom>', ...
    'Color', 'w', ...
    'Name', sprintf('%s motor map', subjectID));

% Render cortical surface.
ieeg_RenderGifti(cortex_R);
ax = gca;
hold(ax, 'on');
alpha(ax, 0.2);

% Electrode coordinates for the cortical-array subset.
els = [locInfoCortec.x(:), ...
       locInfoCortec.y(:), ...
       locInfoCortec.z(:)];

electrodeNames = locInfoCortec.name(:);

% Retrieve signed motor-mapping values for the selected electrodes.
assert(isfield(rvals, r2Field), ...
    'Field "%s" was not found in rvals.', r2Field);
r2 = rvals.(r2Field)(cortecDPIdx);
r2 = r2(:);

assert(size(els, 1) == numel(r2), ...
    'Number of electrodes and r^2 values do not match.');

% Use a symmetric color range centered at zero.
absMaxR2 = max(abs(r2), [], 'omitnan');
if isempty(absMaxR2) || isnan(absMaxR2) || absMaxR2 == 0
    absMaxR2 = 1;
end

r2Plot = max(min(r2, absMaxR2), -absMaxR2);

% Scale marker diameter by the magnitude of r^2. scatter3 interprets marker
% sizes as areas, therefore the visually scaled diameter is squared.
r2SizeScaled = abs(r2Plot) ./ absMaxR2;
elSize = minElSize + (maxElSize - minElSize) .* r2SizeScaled;
elSizeScatter = elSize.^2;

% Exclude electrodes without valid coordinates or r^2 values.
validIdx = ~isnan(r2Plot) & all(~isnan(els), 2);

% Diverging colormap for signed values.
cm = slanCM('Vik');
colormap(ax, cm);
clim(ax, [-absMaxR2 absMaxR2]);

% Black outline layer.
scatter3(ax, ...
    els(validIdx,1), els(validIdx,2), els(validIdx,3), ...
    (elSize(validIdx) + 5).^2, ...
    'k', ...
    'filled');

% Colored electrode layer.
scatter3(ax, ...
    els(validIdx,1), els(validIdx,2), els(validIdx,3), ...
    elSizeScatter(validIdx), ...
    r2Plot(validIdx), ...
    'filled', ...
    'MarkerEdgeColor', 'none');

% Plot all remaining electrodes in black.
ieeg_elAdd([locInfoRest.x, locInfoRest.y, locInfoRest.z], 'k', 15);

% Colorbar.
cb = colorbar(ax);
cb.Label.String = 'signed r^2';
cb.FontSize = 14;
cb.Ticks = [-absMaxR2 0 absMaxR2];
cb.TickLabels = { ...
    sprintf('%.2f', -absMaxR2), ...
    '0', ...
    sprintf('%.2f', absMaxR2)};

% Optional electrode labels. Modify labelNames as needed.
labelNames = {'RQA7','RQA8'};
labelIdx = ismember(strtrim(electrodeNames), labelNames); %#ok<NASGU>

% ieeg_elAddText( ...
%     els(labelIdx,:), ...
%     20, ...
%     15, ...
%     electrodeNames(labelIdx));

axis(ax, 'equal');
axis(ax, 'off');
ieeg_viewLight(0, 90);

%% Render motor map on T1-weighted MRI slices
% Weight each projected electrode by its signed r^2 value and retain only
% statistically significant electrodes.
assert(isfield(rvals, 'p_foot_HFB'), ...
    'Field "p_foot_HFB" was not found in rvals.');
assert(numel(rvals.(r2Field)) == numel(rvals.p_foot_HFB), ...
    'r^2 and p-value vectors must have the same length.');

weights = rvals.(r2Field) .* (rvals.p_foot_HFB < pThresh);

[xSlice, ySlice, zSlice] = seegview_sliceplot( ...
    dpLocs, bvol, x, y, z, sliceThickness);

sv_weight_add(dpLocs, weights, xSlice, threshFrac);
set(gcf, 'Name', sprintf('%s foot %s, sagittal, slice thickness %d mm', ...
    subjectID, taskName, sliceThickness));

sv_weight_add(dpLocs, weights, ySlice, threshFrac);
set(gcf, 'Name', sprintf('%s foot %s, coronal, slice thickness %d mm', ...
    subjectID, taskName, sliceThickness));

sv_weight_add(dpLocs, weights, zSlice, threshFrac);
set(gcf, 'Name', sprintf('%s foot %s, axial, slice thickness %d mm', ...
    subjectID, taskName, sliceThickness));

%% Optionally save MRI-slice figures
if saveSliceFigures
    outputDir = fullfile(dataPath, 'figures', subjectID);
    if ~isfolder(outputDir)
        mkdir(outputDir);
    end

    kjm_printfig(fullfile(outputDir, sprintf('%s_foot_%s_sag_%dmm', ...
        subjectID, taskName, sliceThickness)), 10*[3 2]);

    kjm_printfig(fullfile(outputDir, sprintf('%s_foot_%s_cor_%dmm', ...
        subjectID, taskName, sliceThickness)), 10*[3 2]);

    kjm_printfig(fullfile(outputDir, sprintf('%s_foot_%s_ax_%dmm', ...
        subjectID, taskName, sliceThickness)), 10*[3 2]);
end
