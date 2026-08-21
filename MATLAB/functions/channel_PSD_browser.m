function channel_PSD_browser(spectra, freqs, channelNames, dBAxis)
    % CHANNEL_PSD_BROWSER Interactive GUI for browsing Power Spectral Density (PSD) of EEG channels
    %
    %   channel_PSD_browser(spectra, freqs, channelNames, dBAxis) creates a GUI that allows users
    %   to visualize the PSD of different EEG channels using arrow keys or a text input field.
    %
    %   Inputs:
    %       Required
    %           spectra (MxN matrix) - Power spectral density values, where M is the number of channels and N is the number of frequency bins.
    %           freqs (vector) - Frequency values corresponding to the columns of 'spectra'.
    %       Optional:
    %           channelNames (cell array of strings) - Names of the channels (default: 'Ch 1' to 'Ch M').
    %           dBAxis (logical) - Whether to display the power in decibels (default: true).
    %
    %   Controls:
    %       - Left/Right arrow keys: Navigate between channels.
    %       - Editable text field: Enter a channel number to jump to it.
    %
    %   Dependencies:
    %       - Requires MATLAB R2016b or later for argument validation and GUI features.
    %
    %   Author:
    %       Frederik Lampert, Mayo Clinic , 2025.
    %--------------------------------------------------------------------------
    
    % Validate function arguments
    arguments
        spectra {mustBeVectorOrMatrixDouble}
        freqs {mustBeVector, mustMatchLargerDimension(spectra, freqs)}
        channelNames (1,:) {mustBeText} = arrayfun(@(ch) sprintf('Ch %d', ch), 1:size(spectra,1), 'UniformOutput', false)
        dBAxis (1,1) {mustBeNumericOrLogical} = true
    end
    
    %----------------------------------------------------------
    % Step 1: Initialize GUI Elements
    %----------------------------------------------------------

     % Define channel index boundaries
    minIdx = 1;
    maxIdx = size(spectra,1); 

    % Initialize the GUI
    fig = figure('Name', 'Channel PSD Browser', ...
        'Units', 'Normalized', ...
        'Position', [0.1, 0.15, 0.8, 0.7], ...
        'Color', [0.8, 0.8, 0.8], ...
        'Pointer', 'arrow', ...
        'KeyPressFcn', @keyPressHandler);

    % Create axes for the plot
    ax = axes('Parent', fig, ...
        'Units', 'Normalized', ...
        'InnerPosition', [0.141    0.120    0.8    0.742], ...
        'OuterPosition', [0.025    0.020    0.95    0.911]);

    % Create an editable text field for channel selection
    chEdit = uicontrol('Style', 'edit', ...
        'Units', 'Normalized', ...
        'Position', [ 0.485    0.9    0.06    0.0580], ...
        'String', '1', 'Callback', @textFieldCallback);
    
    % Display usage instructions
    annotation(fig,'textbox',...
    [0.24    0.029    0.516    0.032],...
    'VerticalAlignment','middle',...
    'String','Use Left and Right Arrow Keys to switch between the channels or enter the chnnel number to the editable text box',...
    'HorizontalAlignment','center',...
    'FontSize',12);
    
    %----------------------------------------------------------
    % Step 2: Store Data in Figure for Easy Access
    %----------------------------------------------------------

    % Store data in figure handle
    set(fig,'UserData',struct('Data',spectra, ...
        'Freqs', freqs, ...
        'ChNames', {channelNames}, ...
        'dB', dBAxis, ...
        'CurrentChannel', 1, ...
        'PlotHandle', []));

    % Draw initial PSD plot
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

     % Update the PSD plot based on the selected channel
    function updatePlot();
        ch = fig.UserData.CurrentChannel;
        if isempty(fig.UserData.PlotHandle)
            if fig.UserData.dB
                fig.UserData.PlotHandle = plot(fig.UserData.Freqs, pow2db(fig.UserData.Data(ch,:)), 'color', [0 0.4470 0.7410], 'LineWidth', 2);
                ylabel('Power (dB/Hz)', 'FontSize', 14);
            else
                fig.UserData.PlotHandle = semilogy(fig.UserData.Freqs, fig.UserData.Data(ch,:), 'color', [0 0.4470 0.7410], 'LineWidth', 2);
                ylabel('Power (\muV^2/Hz)', 'FontSize', 14);
            end
            xlabel('Frequencies (Hz)', 'FontSize', 14);
        else
            if fig.UserData.dB
                set(fig.UserData.PlotHandle, 'YData', pow2db(fig.UserData.Data(ch,:)));
            else
                set(fig.UserData.PlotHandle, 'YData', fig.UserData.Data(ch,:));
            end
        end
        grid on
        title(fig.UserData.ChNames{ch}, 'FontSize', 14);
        xlim([fig.UserData.Freqs(1) fig.UserData.Freqs(end)]);
        drawnow;
    end

end

%% Helper Functions
function mustBeVectorOrMatrixDouble(x)
    if ~(isnumeric(x) && ismatrix(x) && isa(x,'double') && ~isempty(x))
        error("'signal' must be a non-empty double vector or matrix.");
    end
end

% Custom validation function to check the size of frequencies
function mustMatchLargerDimension(signal, frequencies)
    signal_size = max(size(signal)); % Get the larger dimension of signal
    if length(frequencies) ~= signal_size
        error("'frequencies' must have the same size as the larger dimension of 'signal'.");
    end
end
