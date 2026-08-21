function chanel_spectogram_browser(signal, sampleRate, window, overlap, frequencies, cmap, flipaxes, channelNames)
%--------------------------------------------------------------------------
% Function: chanel_spectogram_browser
%--------------------------------------------------------------------------
% Description:
%   This function visualizes spectrograms for individual channels of a 
%   multi-channel signal. It provides an interactive GUI to browse through 
%   channels using arrow keys or a text field input.
%
% Usage:
%   chanel_spectogram_browser(signal, sampleRate, window, overlap, 
%                             frequencies, cmap, flipaxes, channelNames)
%
% Input Arguments:
%   Required:
%       - signal (double matrix)         : Multi-channel signal (time x channels)
%       - sampleRate (numeric)           : Sampling frequency in Hz
%   Optional:
%       - window (numeric, optional)     : Window size for spectrogram (default: 5 * sampleRate)
%       - overlap (numeric, optional)    : Overlap between windows (default: floor(window/2))
%       - frequencies (vector, optional) : Frequencies to compute (default: [])
%       - cmap (string/matrix, optional) : Colormap for the spectrogram (default: 'turbo')
%       - flipaxes (logical, optional)   : Whether to flip y-axis (default: true)
%       - channelNames (cell array)      : Names of channels (default: 'Ch 1', 'Ch 2', etc.)
%
% Functionality:
%   - Computes spectrograms for each channel using the specified parameters.
%   - Displays an interactive GUI where users can navigate between channels 
%     using left/right arrow keys or a numeric text field.
%   - Provides a color-mapped visualization of power spectral density in dB.
%
% Dependencies:
%   - Requires Signal Processing Toolbox (for `spectrogram` function).
%
%   Author:
%       Frederik Lampert, Mayo Clinic , 2025.
%--------------------------------------------------------------------------

    arguments
        signal {mustBeVectorOrMatrixDouble}
        sampleRate (1,1) {mustBeNumeric}
        window (1,1) {mustBeNumeric} = 5*sampleRate 
        overlap (1,1) {mustBeNumeric} = floor(window/2)
        frequencies {mustBeVector} = 1:sampleRate/2  
        cmap {mustBeStringOrMatrix} ='turbo'
        flipaxes (1,1) {mustBeNumericOrLogical} = true
        channelNames (1,:) {mustBeText} = arrayfun(@(ch) sprintf('Ch %d', ch), 1:size(signal,2), 'UniformOutput', false)
    end
    
    %----------------------------------------------------------
    % Step 1: Compute Spectrograms
    %----------------------------------------------------------
    disp('Calculating spectograms')
    % spctgrms = cell(1,size(signal,2));

    total_samples = max(size(signal));
    hopSize = window - overlap;
    nCol = fix((total_samples-overlap)/hopSize);                   % Get number of columns
    spctgrms = zeros(length(frequencies),nCol,size(signal,2));     % Initiate empty matrix for spectograms
    for j = 1:size(signal,2)
        [~, f, t,spctgrms(:,:,j)] = spectrogram(signal(:,j),window,overlap,frequencies,sampleRate);
    end

    %----------------------------------------------------------
    % Step 2: Initialize GUI Elements
    %----------------------------------------------------------
    % Specify the lowest and higest possible index
    minIdx = 1;
    maxIdx = size(signal,2); 

    % Initialize the GUI
    fig = figure('Name', 'Channel Spectogrm Browser', ...
        'Units', 'Normalized', ...
        'Position', [0.1, 0.15, 0.8, 0.7], ...
        'Color', [0.8, 0.8, 0.8], ...
        'Pointer', 'arrow', ...
        'KeyPressFcn', @keyPressHandler);

    % Create axes in the middle
    ax = axes('Parent', fig, ...
        'Units', 'Normalized', ...
        'InnerPosition', [0.141    0.120    0.8    0.742], ...
        'OuterPosition', [0.025    0.020    0.95    0.911]);

    % Create editable text field above the axes
    chEdit = uicontrol('Style', 'edit', ...
        'Units', 'Normalized', ...
        'Position', [ 0.47    0.9    0.06    0.0580], ...
        'String', '1', 'Callback', @textFieldCallback);
    
    % Create a texbox with explantionof control
    annotation(fig,'textbox',...
        [0.15    0.02    0.7    0.03],...
        'VerticalAlignment','middle',...
        'String','Use Left and Right Arrow Keys to switch between the channels or enter the chnnel number to the editable text box',...
        'HorizontalAlignment','center',...
        'FontSize',12);
    
    %----------------------------------------------------------
    % Step 3: Store Data in Figure for Easy Access
    %----------------------------------------------------------

    % Save spctogram into figure handles
    set(fig,'UserData', struct('Spectograms', spctgrms, ... 
        'Frequencies', f, ...
        'Time', t, ...
        'CurrentChannel',1, ...
        'ChNames', {channelNames}, ...
        'PlotHandle', []));

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

    % Update the point on the axes
    function updatePlot()
        data = get(fig,'UserData');
        ch = data.CurrentChannel;
        tr = data.Time;
        fr = data.Frequencies;
        pxx = data.Spectograms(:,:,ch);

        if isempty(fig.UserData.PlotHandle)
            fig.UserData.PlotHandle = imagesc(tr,fr,pow2db(abs(pxx)));
            if flipaxes
               axis xy
            end
        else
            set(fig.UserData.PlotHandle, 'CData', pow2db(abs(pxx)));
        end
        xlabel("Time (s)"); ylabel("Frequency (Hz)")
        title(fig.UserData.ChNames{ch});
        colormap(cmap), 
        cb = colorbar;
        cb.Label.String = 'Power (dB/Hz)';
        fontsize(16,'points');
        drawnow;
    end
   
end

%% Helper Functions
function mustBeVectorOrMatrixDouble(x)
    if ~(isnumeric(x) && ismatrix(x) && isa(x,'double') && ~isempty(x))
        error("'signal' must be a non-empty double vector or matrix.");
    end
end

function mustBeStringOrMatrix(x)
    % Accepts either a string or a vector of integers
    % If it's char, convert to string for consistent checks
    if ischar(x)
        x = string(x);
    end

    if isstring(x)
        % valid
        return;
    elseif ismatrix(x)
        % must be a matrix containing RGB values
        if ~any(size(x) == 3)
            error("Colormap must be a N x 3 matrix specifying RGB values ");
        end
    else
        error("specification must be either:\n%s\n%s", ...
            "-> string specifying colormap (""turbo, jet, etc."") or array of RGB values");
    end
end

% Custom validation function to check the size of frequencies
function mustMatchLargerDimension(signal, frequencies)
    signal_size = max(size(signal)); % Get the larger dimension of signal
    if length(frequencies) ~= signal_size
        error("'frequencies' must have the same size as the larger dimension of 'signal'.");
    end
end
