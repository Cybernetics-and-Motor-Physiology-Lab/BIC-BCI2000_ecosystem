function label_add_popout(locs, bad_channels, plot_offset, msize, dispLegend)
% LABEL_ADD_POPOUT Plots electrode positions with a geometry-based popout offset.
%
% ## Description
% This function visualizes electrode positions in 3D space, differentiating
% between good and bad channels. Instead of requiring manually specified
% angular parameters, it computes a local offset direction automatically
% from the electrode trajectory. This is more robust for electrodes with a
% slight curve.
%
% ## Syntax
% ```matlab
% label_add_popout(locs, bad_channels)
% label_add_popout(locs, bad_channels, plot_offset)
% label_add_popout(locs, bad_channels, plot_offset, msize)
% label_add_popout(locs, bad_channels, plot_offset, msize, dispLegend)
% ```
%
% ## Inputs
% - **locs** (`Nx3` matrix): Coordinates of electrode positions. Each row
%   represents an electrode with `[X, Y, Z]` coordinates.
%
% - **bad_channels** (`vector`): Indices of electrodes that are considered
%   "bad" and will be highlighted differently in the plot.
%
% - **plot_offset** (`optional`, `double`): Relative offset applied to
%   electrode positions for visualization purposes. *Default:* `0.08`.
%
% - **msize** (`optional`, `double`): Size of the markers representing electrodes.
%   *Default:* `18`.
%
% - **dispLegend** (`optional`, `logical`): Whether to display the legend.
%   *Default:* `false`.
%
% ## Outputs
% - **None**: The function generates a 3D plot displaying electrode positions
%   with labels. It differentiates between good and bad channels visually.
%
% ## Notes
% - The offset is computed automatically from the local trajectory geometry.
% - For each contact, a local tangent vector is estimated.
% - A perpendicular direction is then computed using a stable reference axis.
% - This works better than manually specifying azimuth/elevation for
%   electrodes with curvature.
%
% -------------------------------------------------------------------------

% ---- Validate inputs ----
if nargin < 2 || isempty(bad_channels)
    bad_channels = [];
end

if nargin < 3 || isempty(plot_offset)
    plot_offset = 0.05;
end

if nargin < 4 || isempty(msize)
    msize = 18;
end

if nargin < 5 || isempty(dispLegend)
    dispLegend = false;
end

if ~isnumeric(locs) || size(locs,2) ~= 3
    error('locs must be an Nx3 numeric matrix.');
end

if ~isnumeric(bad_channels)
    error('bad_channels must be a numeric vector of indices.');
end

bad_channels = bad_channels(:).';
bad_channels = bad_channels(bad_channels >= 1 & bad_channels <= size(locs,1));
bad_channels = unique(round(bad_channels));

%% ---- Compute simple Z-direction offset ----
% Scale offset to the spatial extent of the electrode cloud
scale = plot_offset * max(range(locs, 1));
if scale == 0
    scale = plot_offset;
end

a_offset = repmat([0 0 scale], size(locs,1), 1);

% Apply offset
locs = locs + a_offset;

%% ---- Identify good channels ----
good_chnls = 1:size(locs,1);
good_chnls(bad_channels) = [];

%% ---- Plot good channels ----
g_plts = plot3(locs(good_chnls,1), locs(good_chnls,2), locs(good_chnls,3), 'o', ...
    'MarkerSize', msize, ...
    'LineWidth', 0.5, ...
    'MarkerEdgeColor', 0.4*[1 1 1], ...
    'MarkerFaceColor', 0.99*[1 1 1]);
hold on

%% ---- Plot bad channels ----
if ~isempty(bad_channels)
    b_plts = plot3(locs(bad_channels,1), locs(bad_channels,2), locs(bad_channels,3), 'o', ...
        'MarkerSize', msize, ...
        'LineWidth', 0.5, ...
        'MarkerEdgeColor', 0.4*[1 1 1], ...
        'MarkerFaceColor', [1 0 0]);
else
    b_plts = gobjects(0);
end

%% ---- Add legend ----
if dispLegend
    if ~isempty(bad_channels)
        legend([g_plts, b_plts], {'Good channels', 'Bad channels'}, 'Location', 'best');
    else
        legend(g_plts, {'Good channels'}, 'Location', 'best');
    end
end

%% ---- Label each electrode ----
for k = 1:size(locs,1)
    text(locs(k,1)*1.01, locs(k,2)*1.01, locs(k,3)*1.01, num2str(k), ...
        'FontSize', 10, ...
        'HorizontalAlignment', 'center', ...
        'VerticalAlignment', 'middle');
end

%% ---- Plot formatting ----
xlabel('X-axis');
ylabel('Y-axis');
zlabel('Z-axis');
title('Electrode Positions');
axis equal
hold off

end