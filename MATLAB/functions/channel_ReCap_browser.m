function channel_ReCap_browser(plotData, freqs)
%--------------------------------------------------------------------------
% Function:
%--------------------------------------------------------------------------
% Description:
%   This function visualizes spectrograms for individual channels of a
%   multi-channel signal. It provides an interactive GUI to browse through
%   channels using arrow keys or a text field input.
%
% Usage:
%   channel_ReCap_browser(plotData)
%
% Input Arguments:
%   Required:
%       - plotData (structure)         : structure containing the data
%       - freqs (double)               : array specifiyng frequencies for
%       plotting
%
% Functionality:
%   - Computes spectrograms for each channel using the specified parameters.
%
%   Author:
%       Frederik Lampert, Mayo Clinic , 2025.
%--------------------------------------------------------------------------
if nargin < 2
    freqs = 0:size(plotData(1).PSD,1)-1;
end

%----------------------------------------------------------
% Step 1: Initialize GUI Elements
%----------------------------------------------------------
% Specify the lowest and higest possible index
minIdx = 1;
maxIdx = 32;
chNames = arrayfun(@(ch) sprintf('Ch %d', ch), 1:maxIdx, 'UniformOutput', false);
dates = [plotData.Date];

% Initialize the GUI
fig = figure('Name', 'Channel Recording Capabilities Browser', ...
    'Units', 'Normalized', ...
    'Position', [0.1, 0.15, 0.8, 0.7], ...
    'Color', [0.8, 0.8, 0.8], ...
    'Pointer', 'arrow', ...
    'KeyPressFcn', @keyPressHandler);

% Create editable text field above the axes
chEdit = uicontrol('Style', 'edit', ...
    'Units', 'Normalized', ...
    'Position', [ 0.47    0.95    0.06    0.03], ...
    'String', '1', 'Callback', @textFieldCallback);

% Create a texbox with explantionof control
annotation(fig,'textbox',...
    [0.15    0.02    0.7    0.03],...
    'VerticalAlignment','middle',...
    'String','Use Left and Right Arrow Keys to switch between the channels or enter the chnnel number to the editable text box',...
    'HorizontalAlignment','center',...
    'FontSize',12);

% Create tiled layout
tlo = tiledlayout(fig, 2, 1, 'TileSpacing', 'compact');
ax1 = nexttile(tlo, 1);  % upper: Impedance & RMS
ax2 = nexttile(tlo, 2);  % lower: PSD image

%----------------------------------------------------------
% Step 2: Store Data in Figure for Easy Access
%----------------------------------------------------------

% Store data in figure
set(fig, 'UserData', struct( ...
    'CurrentChannel', 1, ...
    'PlotData', plotData, ...
    'ChannelNames', {chNames}, ...
    'Freqs', freqs, ...
    'Dates', dates, ...
    'Axes', struct('Top', ax1, 'Bottom', ax2)));

% Draw initial PSD
updatePlot();

%----------------------------------------------------------
% Callback Functions
%----------------------------------------------------------

% Callback function for keypress events
    function keyPressHandler(~, event)
        switch event.Key
            case 'leftarrow'
                fig.UserData.CurrentChannel = max(minIdx, fig.UserData.CurrentChannel - 1); % Decrement x, ensure it doesn't go below 1
            case 'rightarrow'
                fig.UserData.CurrentChannel = min(maxIdx, fig.UserData.CurrentChannel + 1); % Increment x, ensure it doesn't exceed32
        end
        updatePlot();
        set(chEdit, 'String', num2str(fig.UserData.CurrentChannel)); % Update the text field
    end

% Callback function for the editable text field
    function textFieldCallback(src, ~)
        inputValue = str2double(get(src, 'String'));
        if isnan(inputValue) || inputValue < minIdx || inputValue > maxIdx
            errordlg(sprintf('Please enter an integer between %d and %d.', minIdx, maxIdx), 'Invalid Input');
            set(src, 'String', num2str(ch)); % Reset to the current value of x
        else
            fig.UserData.CurrentChannel = round(inputValue); % Update x with the new value
            updatePlot();
        end
    end

% --- Plotting logic
    function updatePlot()
        ch = fig.UserData.CurrentChannel;
        pd = fig.UserData.PlotData;
        dates = fig.UserData.Dates;
        nDays = numel(dates);

        % Clear previous content
        cla(fig.UserData.Axes.Top);
        cla(fig.UserData.Axes.Bottom);

            % -- Sublot 1: Dual y-axis plot for Impedance and RMS
            yyaxis(fig.UserData.Axes.Top, 'left');
            imp = cell2mat(arrayfun(@(x) x.Impedance(ch), pd, 'UniformOutput', false));
            plot(fig.UserData.Axes.Top, dates, imp, '-o', 'Color', [0.2 0.4 0.8], 'LineWidth', 1.5);
            ylabel(fig.UserData.Axes.Top,'Impedance [Ohms]');
    
            yyaxis(fig.UserData.Axes.Top, 'right');
            rms = cell2mat(arrayfun(@(x) x.RMS(ch), pd, 'UniformOutput', false));
            plot(fig.UserData.Axes.Top, dates, rms, '-s', 'Color', [0.8 0.2 0.4], 'LineWidth', 1.5);
            ylabel(fig.UserData.Axes.Top,'PSD RMS');
    
            title(fig.UserData.Axes.Top, ['Channel ', num2str(ch), ' - Impedance & RMS']);
            xlim(fig.UserData.Axes.Top,'tight')
            xlabel('Date');
            grid(fig.UserData.Axes.Top, 'on');

            % -- Sublot 2: spectogram
            PSDs = cat(2, pd.PSD);
            PSDs = PSDs(:, ch:maxIdx:end);  % extract ch from all days
            imagesc(fig.UserData.Axes.Bottom, dates, fig.UserData.Freqs, pow2db(PSDs));
            set(fig.UserData.Axes.Bottom, 'YDir', 'normal');
            ylabel(fig.UserData.Axes.Bottom,'Frequency (Hz)');
            xlabel('Date');
            cb = colorbar(fig.UserData.Axes.Bottom);
            cb.Label.String = 'Power (dB/Hz)';
            title(fig.UserData.Axes.Bottom, sprintf('Channel %d - PSD Evolution', ch));
            colormap(fig.UserData.Axes.Bottom, turbo);
            fontsize(16,"points")
    end
end

