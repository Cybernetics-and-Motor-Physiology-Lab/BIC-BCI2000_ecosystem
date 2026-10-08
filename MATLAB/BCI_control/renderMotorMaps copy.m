clear all; close all; clc
% Render brain for patient IMG
% Jordan Bilderbeek, Frederik Lampert June 5, 2026

%% Add Paths
addpath('./functions/');
addpath(genpath('.../mnl_seegview')); % https://github.com/MultimodalNeuroimagingLab/mnl_seegview
addpath(genpath('path/to/mnl_ieegBasics/functions')); % https://github.com/MultimodalNeuroimagingLab/mnl_ieegBasics
addpath(genpath('.../slanCM')); % https://www.mathworks.com/matlabcentral/fileexchange/120088-200-colormaps 

%% Load /specify paths to data
% Electrode locations and brain render
electrodes_tsv_name = '../data/BCI/IMG_electrodes.tsv';

% Load r2 values from motor mapping
load('../data/BCI/IMG_r_mot_rvals_noBon.mat');
load('../data/brainMAT/sub-IMG_brain.mat')

%% load electrodes + surface + meshes
loc_info = ieeg_readtableRmHyphens(electrodes_tsv_name);

locs = [loc_info.x, loc_info.y, loc_info.z];
[dp_channels, dp_locs] = locs_DPRR(locs,5,'y');

loc_info_DP = loc_info(dp_channels, :);

%% RSI10-11 is 34-35 | for d3s2 RQA7-8 is 120-121

loc_info_cortec = loc_info_DP; 
cortec_DP_idx = [74:101]; % d3s2 RQA7-8 is 120-121 | d4s2 (imagery BCI) RSI10-11 - 27:54
loc_info_cortec = loc_info_cortec(cortec_DP_idx, :); 

loc_info_rest =  loc_info_DP;
loc_info_rest(cortec_DP_idx, :) = []; 

gR       = gifti(gii_name);

%% figure

fig = figure('Position',[100 100 2000 1500], 'Visible','on', ...
    'PaperUnits','inches', 'PaperPosition',[0,0,20,15], ...
    'PaperPositionMode','auto', 'PaperSize',[20,15], ...
    'PaperType','<custom>');

% Render cortex
ieeg_RenderGifti(gR);
ax = gca;
hold(ax, 'on');
alpha(.2);

% Electrode coordinates
els = [loc_info_cortec.x(:), ...
       loc_info_cortec.y(:), ...
       loc_info_cortec.z(:)];

electrode_names = loc_info_cortec.name(:);

% Your signed values
r2 = rvals.r_foot_HFB(cortec_DP_idx);%r2 = rvals.r_hand_HFB(27:54);
r2 = r2(:);

% Safety check
assert(size(els,1) == numel(r2), 'Number of electrodes and r2 values do not match.');

min_elsize = 5;
max_elsize = 20;

absmax_r2 = max(abs(r2), [], 'omitnan');

if isempty(absmax_r2) || isnan(absmax_r2) || absmax_r2 == 0
    absmax_r2 = 1;
end

% Clip values to symmetric range
r2_plot = max(min(r2, absmax_r2), -absmax_r2);

% Size is based on magnitude only
r2_size_scaled = abs(r2_plot) ./ absmax_r2;
el_size = min_elsize + (max_elsize - min_elsize) .* r2_size_scaled;

% scatter3 sizes are area-like, so square them for stronger visual scaling
el_size_scatter = el_size.^2;

% Valid electrodes only
valid_idx = ~isnan(r2_plot) & all(~isnan(els), 2);

% slan CM as dependency 
cm = slanCM('Vik');

colormap(ax, cm);
caxis(ax, [-absmax_r2 absmax_r2]);

% Black outline layer
scatter3(ax, ...
    els(valid_idx,1), els(valid_idx,2), els(valid_idx,3), ...
    (el_size(valid_idx) + 5).^2, ...
    'k', ...
    'filled');

% Colored electrode layer
scatter3(ax, ...
    els(valid_idx,1), els(valid_idx,2), els(valid_idx,3), ...
    el_size_scatter(valid_idx), ...
    r2_plot(valid_idx), ...
    'filled', ...
    'MarkerEdgeColor', 'none');

ieeg_elAdd([loc_info_rest.x, loc_info_rest.y, loc_info_rest.z], 'k', 15)


cb = colorbar(ax);
cb.Label.String = 'signed r^2';
cb.FontSize = 14;
cb.Ticks = [-absmax_r2 0 absmax_r2];

% Optional: make the tick labels explicit
cb.TickLabels = { ...
    sprintf('%.2f', -absmax_r2), ...
    '0', ...
    sprintf('%.2f', absmax_r2)};


label_idx = ismember(strtrim(electrode_names), {'RQA7','RQA8'});

% ieeg_elAddText( ...
%     els(label_idx,:), ...
%     20, ...
%     15, ...
%     electrode_names(label_idx));


axis(ax, 'equal');
axis(ax, 'off');
set(ax, 'CLim', [-absmax_r2 absmax_r2]);

ieeg_viewLight(0,90)

%% SEEG view  plot on T1 MRI slices
pt = 'IMG';
pt_task = 'overtBCI_d3s2';


   side = 'B';
   figsave = 'n';
   slicethickness = 5;

    % save figs?
    figsave = 0;
    
    % set parameters for slice plots and rendering angles
    slicethickness = 5; % thickness of each plotted slice, in mm
    mr_clims=[0 1]; % adjusts brightness of plots
    pthresh= 0.05; % defines which rvals are significant
    threshFrac = 0.01; % abs(values) below threshFrac*wmaximum weight are not plotted with color
    
    % Plot slice of brain image with weighted, unlabelled electrodes
    wts = rvals.r_foot_HFB.*(rvals.p_foot_HFB<pthresh); % set weights F
    [x_slice, y_slice, z_slice] = seegview_sliceplot(dp_locs, bvol, x, y, z, slicethickness,mr_clims, side);
    sv_weight_add(dp_locs, wts, x_slice, threshFrac); set(gcf,'Name',[pt ' foot ' pt_task ', sagittal, slice thickness ' num2str(slicethickness) 'mm']);
    sv_weight_add(dp_locs, wts, y_slice, threshFrac); set(gcf,'Name',[pt ' foot ' pt_task ', coronal, slice thickness ' num2str(slicethickness) 'mm']);
    sv_weight_add(dp_locs, wts, z_slice, threshFrac); set(gcf,'Name',[pt ' foot ' pt_task ', axial, slice thickness ' num2str(slicethickness) 'mm']);
    
    if figsave
        axes(x_slice.ha(1)); kjm_printfig(['data/' pt '/' pt '_foot_2_sec_nobon_' pt_task '_sag_' num2str(slicethickness) 'mm'],10*[3 2])
        axes(y_slice.ha(1)); kjm_printfig(['data/' pt '/' pt '_foot_2_sec_nobon_' pt_task '_cor_' num2str(slicethickness) 'mm'],10*[3 2])
        axes(z_slice.ha(1)); kjm_printfig(['data/' pt '/' pt '_foot_2_sec_nobon_' pt_task '_ax_' num2str(slicethickness) 'mm'],10*[3 2])
    end