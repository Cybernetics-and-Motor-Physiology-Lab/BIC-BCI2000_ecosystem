function varargout = channel_trialspect_browser(spectra, channelNames, colorlims, clrbar_title)
    arguments
        spectra (:,:,:) {mustBeNumeric}
        channelNames cell = {} % optional: cell array of strings/chars
        colorlims  double = [] % optional: vector [min max]
        clrbar_title string = "" % optional: string title
    end   

    % Default assignment if channelNames is empty
    if isempty(channelNames)
        channelNames = arrayfun(@(ch) sprintf('Ch %d', ch), 1:size(spectra,1), 'UniformOutput', false);
    end

    % Specify the lowest and higest possible index
    minIdx = 1;
    maxIdx = size(spectra,1); 

    % Initialize the GUI
    fig = figure('Name', 'Channel Browser', ...
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
        'Position', [ 0.5410    0.8680    0.0590    0.0580], ...
        'String', '1', 'Callback', @textFieldCallback);
    
    % Create a texbox with explantionof control
    annotation(fig,'textbox',...
    [0.24    0.029    0.516    0.032],...
    'VerticalAlignment','middle',...
    'String','Use Left and Right Arrow Keys to switch between the channels or enter the chnnel number to the editable text box',...
    'HorizontalAlignment','center',...
    'FontSize',12);
   
    % Initialize default channel to dispay (deafault ch=1)
    ch = 1;

    % Draw initial trial spectra
    Himg = imagesc(ax, squeeze(spectra(ch,:,:)));
    title(channelNames{ch}, 'FontSize',14);
    xlabel('#Trial', 'FontSize',14); ylabel('Frequency [Hz]', 'FontSize',14);
    axis xy

    % Create colorbar
    c = colorbar(ax);

    if isempty(colorlims)
        clim(ax,'auto')
    else 
        clim(ax,colorlims)
    end

    if ~isempty(clrbar_title)
        c.Label.String = clrbar_title; 
        c.Label.FontSize = 14;
    end


    % Callback function for keypress events
    function keyPressHandler(~, event)
        switch event.Key
            case 'leftarrow'
                ch = max(minIdx, ch - 1); % Decrement x, ensure it doesn't go below 1
            case 'rightarrow'
                ch = min(maxIdx, ch + 1); % Increment x, ensure it doesn't exceed32
        end
        updatePlot();
        set(chEdit, 'String', num2str(ch)); % Update the text field
    end

    % Callback function for the editable text field
    function textFieldCallback(src, ~)
        inputValue = str2double(get(src, 'String'));
        if isnan(inputValue) || inputValue < minIdx || inputValue > maxIdx
            errordlg(sprintf('Please enter an integer between %d and %d.', minIdx, maxIdx), 'Invalid Input');
            set(src, 'String', num2str(ch)); % Reset to the current value of x
        else
            ch = round(inputValue); % Update x with the new value
            updatePlot();
        end
    end

    % Update the point on the axes
    function updatePlot();
        set(Himg, 'CData', squeeze(spectra(ch,:,:))); % Update the data in the plot
        title(channelNames{ch}, 'FontSize',14);
        drawnow;
    end


end