function varargout = segment_signal(signal, cue, srate, fname)

% segment_signal(signal, cue, srate)
%
% This function creates the GUI that can be used to segment your input
% signal timeseries into different epochs.
%
% This can be EMG, vibrotactile sensor channels etc.
%
% Input parameters:
%   a) signal : time x channels
%   b) cue    : screen cue vector (optional)
%   c) srate  : sampling rate
%
% Click boxes:
%   a) 'Add segment'    : will create a new row of two editable boxes
%       - left box      : insert integer values corresponding to the
%                         marking line you are using, ie 1, 2, 3 etc.
%       - right box     : type the modality or body part you are marking, ie
%                         Hand, Tongue, Foot etc.
%   b) 'Delete segment' : will delete the last row of editable boxes.
%
% Sliders:
%   a) right  : zoom into trace.
%   b) bottom : move trace forward when zoomed in
%
% Marking lines:
%   Use number buttons. Every number has its own color.
%
% Output parameters. Click on 'Save data'.
%   Please note that a warning will appear if the number of marking
%   lines is not even.
%
%   a) 'Modalities'   : The content of the editable boxes will be saved
%                       in a cell format for future reference.
%   b) 'Marked ids'   : Vector with the sequence of marking lines.
%   c) 'Marked inds'  : Vector corresponding to the time points of marked
%                       lines.
%   d) 'Beh':         : Behavioral vector
%
%   It will also save the original input signal.
%
% Panos Kerezoudis, CaMP lab, 2024.
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% INITIALIZING WAIT BOX AND FIGURE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Assign an empty array or default value if cue is not provided
if nargin < 2 || isempty(cue), cue = []; end

% uiwait(msgbox([{'This is the GUI that allows you to epoch your timeseries data.', 'Let''s get started!'}], ...
%         'Getting Started'));

% Initialize the main figure window
% Initialize the main figure window
MAIN = figure('Units', 'Normalized', 'Position', [0.03, 0.05, 0.95, 0.8], ...
    'Color', [0.8, 0.8, 0.8], 'MenuBar', 'none', 'Pointer', 'arrow', ...
    'Visible', 'on', 'Name', 'Signal Viewer', 'NumberTitle', 'off', ...
    'WindowButtonMotionFcn', @mouseHoverCallback, ...
    'WindowKeyPressFcn', @keyPressCallback, ...
    'WindowButtonDownFcn', @mouseClickCallback);

% Initialize the GUI and store the signal, originalName, and marked_inds and marked_ids
handles.signal = signal;
handles.originalName = inputname(1);  % Store the original variable name
handles.marked_inds = [];
handles.marked_ids = [];
handles.modalities = {};
handles.xline_handles = [];
handles.srate = srate;

% Save the handles structure
guidata(MAIN, handles);

% Time vector based on the signal and sample rate
t = (1:length(signal)) / srate;

% Create an axes for the large plot on the left
signal_axes = axes('Parent', MAIN, 'Units', 'Normalized', ...
    'Position', [0.025, 0.1, 0.8, 0.85], 'XColor', 'k', 'YColor', 'k', 'Box', 'on');

% Adjust xlim to match the signal's time range
xlim(signal_axes, [t(1), t(end)]);

% Check if file name was specified
if nargin > 3 && exist('fname','var')
    fig_title = fname;
    strrep(fig_title,'_',' ');
else
    fig_title = 'Input Signal and Cue';
end

% Plot the multichannel signal on a single set of axes with normalization
plot_multichannel_signal(signal, t, signal_axes, cue, fig_title);

% Initial zoom and pan parameters
zoom_factor = 1; % Default zoom
pan_limit = 0;

% Create crosshair lines (hidden by default)
vertical_line = line(signal_axes, [0, 0], ylim(signal_axes), 'Color', 'k', ...
    'LineWidth', 1, 'Visible', 'off');
horizontal_line = line(signal_axes, xlim(signal_axes), [0, 0], 'Color', 'k', ...
    'LineWidth', 1, 'Visible', 'off');

% Colors for xlines based on the integer (1 to 9)
xline_colors = {[0, 0.3, 0], [0, 0, 1], [1, 0, 1], [0.6, 0.3, 0], [0.6, 1, 0.6], ...
    [0.7, 0.9, 1], [1, 0.5, 0], [0, 0.5, 0.5], [0.5 0.5 0.5]};

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% GUI boxes %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Define initial position for editable boxes and pushbuttons
editbox_width1 = 0.02;   % Width of left editable box
editbox_width2 = 0.07;   % Width of right editable box
editbox_height = 0.05;   % Height of editable boxes
box_spacing = 0.02;      % Spacing between boxes
initial_box_y = 0.75;    % Initial Y position for the first row of boxes

% Add "Add Box" pushbutton using the helper function
addBoxButton = uicontrol('Style', 'pushbutton', 'Units', 'Normalized', ...
    'Position', [0.88, initial_box_y + editbox_height + box_spacing, 0.10, editbox_height], ...
    'BackgroundColor', [0.678, 0.847, 0.902], 'String', 'Add Segment', ...
    'FontSize', 14, 'FontWeight', 'bold', 'Callback', @addNewBoxes);

% Add "Delete Box" pushbutton to delete the last added editable box
deleteBoxButton = uicontrol('Style', 'pushbutton', 'Units', 'Normalized', ...
    'Position', [0.88, initial_box_y + (2 * editbox_height) + (2 * box_spacing), 0.10, editbox_height], ...
    'BackgroundColor', [1, 0.6, 0.6], 'String', 'Delete Segment', ...
    'FontSize', 14, 'FontWeight', 'bold', 'Callback', @deleteLastBox);

% Add a yellow "Save Data" button using the helper function
saveButton = uicontrol('Style', 'pushbutton', 'Units', 'Normalized', ...
    'Position', [0.88, 0.13, 0.10, 0.05], ...
    'BackgroundColor', [1, 1, 0.6], 'String', 'Save Data', ...
    'FontSize', 14, 'FontWeight', 'bold', 'Callback', @saveDataCallback);

% Add a Green "Export & Continue" button using the helper function
exportButton = uicontrol('Style', 'pushbutton', 'Units', 'Normalized', ...
    'Position', [0.88, 0.05, 0.10, 0.05], ...
    'BackgroundColor', [0.39 0.83 0.07], 'String', 'Export & Continue', ...
    'FontSize', 14, 'FontWeight', 'bold', 'Callback', @exportDataCallback);

% Variable to store the number of additional rows of boxes
num_boxes = 0;

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Slider bars %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function sliderHandle = createSlider(pos, callback)
        sliderHandle = uicontrol('Parent', MAIN, 'Style', 'slider', 'Units', ...
            'normalized', 'Position', pos, 'Callback', callback); end

% Vertical Zoom slider
zoomSlider = createSlider([0.82, 0.1, 0.03, 0.80], @zoomSliderCallback);

ZoomText = annotation('textbox', [0.825, 0.88, 0.05, 0.05], ...
    'Units', 'normalized', 'Color', 'black', 'FontSize', 14, ...
    'EdgeColor', 'none', 'FontWeight', 'bold', ...
    'String', 'Zoom (>1)');

% Horizontal pan slider
panSlider = createSlider([0.025, 0.02, 0.8, 0.03], @panSliderCallback);

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Pushbuttons %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Function to add new editable boxes
    function addNewBoxes(~, ~)
        % Get the handles structure to access stored variables
        handles = guidata(gcf);

        % Update the number of boxes
        if ~isfield(handles, 'num_boxes')
            handles.num_boxes = 0;  % Initialize num_boxes if it doesn't exist
        end
        handles.num_boxes = handles.num_boxes + 1;

        % Calculate the new Y position for the next row of editable boxes
        new_box_y = initial_box_y - ((handles.num_boxes - 1) * (editbox_height + box_spacing));

        % Create new editable boxes in the next row
        % Left box with integer validation
        left_box = uicontrol('Style', 'edit', 'Units', 'Normalized', ...
            'Position', [0.88, new_box_y, editbox_width1, editbox_height], ...
            'Parent', MAIN, 'Callback', @validateInteger, 'Tag', num2str(handles.num_boxes * 2 + 1), ...
            'FontSize', 12, 'FontWeight', 'bold');

        % Right box without integer validation
        right_box = uicontrol('Style', 'edit', 'Units', 'Normalized', ...
            'Position', [0.91, new_box_y, editbox_width2, editbox_height], ...
            'Parent', MAIN, 'Tag', num2str(handles.num_boxes * 2 + 2), ...
            'FontSize', 12, 'FontWeight', 'bold');

        % Store the handles in the cell array, each row as {left_box, right_box}
        if ~isfield(handles, 'editableBoxHandles')
            handles.editableBoxHandles = {};  % Initialize editableBoxHandles if it doesn't exist
        end
        handles.editableBoxHandles{handles.num_boxes, 1} = left_box;
        handles.editableBoxHandles{handles.num_boxes, 2} = right_box;

        % Update the handles structure
        guidata(gcf, handles);
    end

% Function to delete the last editable boxes
    function deleteLastBox(~, ~)
        % Get the handles structure to access stored variables
        handles = guidata(gcf);

        if handles.num_boxes > 0
            % Find the tags of the last two editable boxes
            edit1_tag = num2str(handles.num_boxes * 2 + 1);
            edit2_tag = num2str(handles.num_boxes * 2 + 2);

            % Delete the last added editable boxes
            delete(findobj('Tag', edit1_tag));
            delete(findobj('Tag', edit2_tag));

            % Remove the last added editable box handles
            handles.editableBoxHandles(handles.num_boxes, :) = [];

            % Decrease the count of boxes
            handles.num_boxes = handles.num_boxes - 1;

            % Recalculate the positions of remaining boxes if necessary
            for i = 1:handles.num_boxes
                new_box_y = initial_box_y - ((i - 1) * (editbox_height + box_spacing));
                set(handles.editableBoxHandles{i, 1}, 'Position', [0.88, new_box_y, editbox_width1, editbox_height]);
                set(handles.editableBoxHandles{i, 2}, 'Position', [0.91, new_box_y, editbox_width2, editbox_height]);
            end
        end

        % Update the handles structure
        guidata(gcf, handles);
    end


% Callback function to validate integer input
    function validateInteger(src, ~)
        input_value = str2double(get(src, 'String'));

        % Check if the input is an integer
        if isnan(input_value) || mod(input_value, 1) ~= 0
            % Display error message if not an integer
            errordlg('Please enter an integer value.', 'Invalid Input', 'modal');
            % Clear the invalid input
            set(src, 'String', '');
        end
    end

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Designate epochs %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function keyPressCallback(src, event)
        % Get the handles structure to access stored variables
        handles = guidata(src);

        % Check if the pressed key corresponds to numbers 1 to 9 or 'x'
        if ismember(event.Key, {'1', '2', '3', '4', '5', '6', '7', '8', '9', 'x'})
            % Get the current axes and mouse position
            ax = gca;  % Get the current axes handle
            mouse_pos = get(ax, 'CurrentPoint');
            x_pos = mouse_pos(1, 1);

            % Ensure the x position is within the valid time range
            if x_pos >= t(1) && x_pos <= t(end)
                % Get the integer value or detect the 'x' key
                if strcmp(event.Key, 'x')
                    % Mark excluded region with red line
                    h = xline(ax, x_pos, 'Color', [1, 0, 0], 'LineWidth', 2);
                    handles.marked_ids(end + 1) = -1;  % Exclude marker
                else
                    % Mark with corresponding color for keys 1 to 9
                    int_val = str2double(event.Key);
                    xline_colors = {[0, 1, 0], [0, 0, 1], [1, 0, 1], [0.6, 0.3, 0], ...
                        [0.6, 1, 0.6], [0.7, 0.9, 1], [1, 0.5, 0], [0, 0.5, 0.5], [0.5, 0.5, 0.5]};
                    color = xline_colors{int_val};
                    h = xline(ax, x_pos, 'Color', color, 'LineWidth', 2);
                    handles.marked_ids(end + 1) = int_val;  % Save the pressed key value
                end
                % Save the marked index and store the xline handle
                handles.marked_inds(end + 1) = x_pos;
                handles.xline_handles(end + 1) = h;  % Store the xline handle

                % Update the handles structure with the new values
                guidata(src, handles);
            end
        end

        % Re-focus the figure to ensure keypress capture
        figure(src);  % Bring the figure back to focus
    end

% Remove nearest xline
    function mouseClickCallback(src, ~)
        % Get the handles structure to access stored variables
        handles = guidata(src);

        % Check if the right mouse button is clicked
        if strcmp(get(src, 'SelectionType'), 'alt')  % 'alt' is the right-click event
            ax = gca;
            mouse_pos = get(ax, 'CurrentPoint');
            x_pos = mouse_pos(1, 1);  % Get the x position of the click

            % Ensure the xline handles are valid and not empty
            if isfield(handles, 'xline_handles') && ~isempty(handles.xline_handles)
                closest_distance = inf;
                closest_idx = -1;

                % Iterate over the stored xline handles
                for i = 1:length(handles.xline_handles)
                    if ishandle(handles.xline_handles(i))  % Check if it's a valid handle
                        xline_pos = get(handles.xline_handles(i), 'Value');  % Get x position of the xline
                        distance = abs(xline_pos - x_pos);  % Compute distance from click to xline
                        if distance < closest_distance
                            closest_distance = distance;
                            closest_idx = i;  % Store the index of the closest xline
                        end
                    end
                end

                % If a close enough xline is found, delete it
                if closest_idx > 0
                    delete(handles.xline_handles(closest_idx));  % Delete the nearest xline
                    handles.xline_handles(closest_idx) = [];  % Remove the handle from the list
                    handles.marked_inds(closest_idx) = [];  % Remove the corresponding marked index
                    handles.marked_ids(closest_idx) = [];  % Remove the corresponding marked id
                end
            end

            % Update handles structure
            guidata(src, handles);
        end
    end


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Fns related to visualization %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Callback function for mouse hover detection with crosshairs
    function mouseHoverCallback(~, ~)
        % Get the current point relative to the axes
        mouse_pos = get(signal_axes, 'CurrentPoint');
        x_pos = mouse_pos(1, 1);
        y_pos = mouse_pos(1, 2);

        % Get current axes limits
        x_limits = xlim(signal_axes);
        y_limits = ylim(signal_axes);

        % Check if the mouse is inside the axes limits
        if x_pos >= x_limits(1) && x_pos <= x_limits(2) && y_pos >= y_limits(1) && y_pos <= y_limits(2)
            set(MAIN, 'Pointer', 'crosshair');  % Change to crosshair

            % Update crosshair lines and make them visible
            set(vertical_line, 'XData', [x_pos, x_pos], 'YData', y_limits, 'Visible', 'on');
            set(horizontal_line, 'XData', x_limits, 'YData', [y_pos, y_pos], 'Visible', 'on');
        else
            set(MAIN, 'Pointer', 'arrow');  % Revert to default arrow

            % Hide the crosshair lines
            set(vertical_line, 'Visible', 'off');
            set(horizontal_line, 'Visible', 'off');
        end
    end

% Zoom Slider Callback: zoom in and out
    function zoomSliderCallback(src, ~)
        zoom_value = get(src, 'Value');
        zoom_factor = 1 + (10 * zoom_value); % Adjust zoom factor as needed
        zoomInOut();
    end

% Pan Slider Callback: move the plot left or right
    function panSliderCallback(src, ~)
        pan_limit = get(src, 'Value');
        zoomInOut();  % Apply panning together with zoom
    end

% Function to zoom in/out and pan the plot
    function zoomInOut()
        total_time = t(end) - t(1);
        zoom_range = total_time / zoom_factor;

        % Calculate new x-limits for zoom
        new_xlim_start = t(1) + pan_limit * (total_time - zoom_range);
        new_xlim_end = new_xlim_start + zoom_range;

        % Ensure the limits stay within the data range
        if new_xlim_start < t(1)
            new_xlim_start = t(1);
            new_xlim_end = new_xlim_start + zoom_range;
        end
        if new_xlim_end > t(end)
            new_xlim_end = t(end);
            new_xlim_start = new_xlim_end - zoom_range;
        end

        % Apply the new x-limits
        xlim(signal_axes, [new_xlim_start, new_xlim_end]);

        % Update crosshair horizontal line limit
        set(horizontal_line, 'XData', xlim(signal_axes));
    end

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Save data %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function saveDataCallback(~, ~)
        % Get the handles structure to access stored variables
        handles = guidata(gcbo);  % Retrieve the latest version of handles

        % Ensure there is an even number of marked_ids to compare
        if mod(length(handles.marked_ids), 2) == 0
            % Compare alternating elements of marked_ids
            if isequal(handles.marked_ids(1:2:end), handles.marked_ids(2:2:end))
                beh = zeros(size(handles.signal, 1), 1);
                for n = 1:2:length(handles.marked_ids)
                    beh(round(handles.marked_inds(n)*handles.srate):round(handles.marked_inds(n+1)*handles.srate)) = handles.marked_ids(n);
                end
            end
        else
            warndlg('The number of marked lines is odd, please doublecheck', 'Warning');
            return;  % Exit the function to prevent saving
        end

        % Retrieve the signal and originalName from the handles
        signal = handles.signal;
        originalName = handles.originalName;

        % Extract marked_ids and marked_inds from handles into local variables
        marked_ids = handles.marked_ids;
        marked_inds = handles.marked_inds;

        % Open file dialog to choose where to save
        [file, path] = uiputfile('*.mat', 'Save Data As');

        if isequal(file, 0)
            disp('User canceled save.');
        else
            % Initialize a cell array to store the box content (renamed to modalities)
            modalities = cell(handles.num_boxes, 2);  % Two columns: left box content, right box content

            % Loop through each pair of editable boxes and get their content
            for i = 1:handles.num_boxes
                left_box_value = get(handles.editableBoxHandles{i, 1}, 'String');  % Get left box content
                right_box_value = get(handles.editableBoxHandles{i, 2}, 'String');  % Get right box content

                % Store the content in the cell array
                modalities{i, 1} = left_box_value;
                modalities{i, 2} = right_box_value;
            end

            % Initialize the structure to store data dynamically
            dataStruct = struct();  % Make sure the structure is initialized

            % Create a field with the original variable name and assign the signal
            dataStruct.(originalName) = signal;

            % Save the signal separately using its original name
            save(fullfile(path, file), '-struct', 'dataStruct', originalName);

            % Save the other variables (marked_ids, marked_inds, and modalities)
            save(fullfile(path, file), 'marked_ids', 'marked_inds', 'srate', 'beh', 'modalities', '-append');

            varargout = {}; % No output needed for this option

            % Display message after successful save
            msgbox('Data has been successfully saved!', 'Success');
            disp(['Data saved to: ', fullfile(path, file)]);
        end
    end

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Export data to workspace %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function exportDataCallback(~, ~)
        % Get the handles structure to access stored variables
        handles = guidata(gcbo);  % Retrieve the latest version of handles

        % Ensure there is an even number of marked_ids to compare
        if mod(length(handles.marked_ids), 2) == 0
            % Compare alternating elements of marked_ids
            if isequal(handles.marked_ids(1:2:end), handles.marked_ids(2:2:end))
                beh = zeros(size(handles.signal, 1), 1);
                for n = 1:2:length(handles.marked_ids)
                    beh(round(handles.marked_inds(n)*handles.srate):round(handles.marked_inds(n+1)*handles.srate)) = handles.marked_ids(n);
                end
            end
        else
            warndlg('The number of marked lines is odd, please doublecheck', 'Warning');
            return;  % Exit the function to prevent saving
        end


        % Extract marked_ids and marked_inds from handles into local variables
        marked_ids = handles.marked_ids;
        marked_inds = handles.marked_inds;

        % Initialize a cell array to store the box content (renamed to modalities)
        modalities = cell(handles.num_boxes, 2);  % Two columns: left box content, right box content

        % Loop through each pair of editable boxes and get their content
        for i = 1:handles.num_boxes
            left_box_value = get(handles.editableBoxHandles{i, 1}, 'String');  % Get left box content
            right_box_value = get(handles.editableBoxHandles{i, 2}, 'String');  % Get right box content

            % Store the content in the cell array
            modalities{i, 1} = left_box_value;
            modalities{i, 2} = right_box_value;
        end
        handles.modalities = modalities;

        % Assign variables to base workspace
        assignin('base', 'marked_inds', marked_inds);
        assignin('base', 'marked_ids', marked_ids);
        assignin('base', 'modalities', modalities);

        % Display message after successful save
        msgbox('Data has been successfully annotated!', 'Success');
        disp('Annotations exported to workspace');
        uiresume(MAIN); % Resume execution to return outputs
    end

% Wait for user interaction and retrieve outputs if Export is clicked
uiwait(MAIN);
if nargout > 0
    varargout = {handles.marked_inds, handles.marked_ids, handles.modalities}; % Assign outputs if needed
else
    varargout = {};
end
close(MAIN)

end


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Embedded plotting fn %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function plot_multichannel_signal(signal, t, signal_axes, cue, fig_title)
    num_chan = size(signal, 2);

    % Normalize each channel independently to scale between a larger range [0, 2]
    for j = 1:num_chan
        signal(:, j) = 3 * (signal(:, j) - min(signal(:, j))) / (max(signal(:, j)) - min(signal(:, j)));
    end

    % Vertical offset between each channel to bring them closer together
    % Adjusting this value will bring the traces closer together vertically
    y_offset = 2; % Adjust to fit better within the plot

    hold(signal_axes, 'on');

    % Plot each channel at a different vertical position, ensuring they stay centered
    for j = 1:num_chan
        plot(signal_axes, t, signal(:, j) + (num_chan - j) * y_offset, 'k', 'LineWidth', 1.5);
    end

    % Get the current y-limits of the axes
    curr_ylim = ylim(signal_axes);

    % Define a color map with slight variations of light orange for the cue
    % Add more colors depending on screen cue
    color_map = [
        1, 0.65, 0;     % Light orange
        1, 0.7, 1;    % Slightly more red-orange
        0.5, 0.5, 0.5;      % Slightly darker orange
        1, 0.75, 0.4    % Lighter orange with more yellow
    ];

    % Find the regions where the cue is non-zero (active)
    if nargin > 3 && ~isempty(cue)  % Check if cue is provided
        cue_active = cue ~= 0;  % Boolean vector, true where cue is non-zero
        transitions = diff([0; cue_active(:); 0]);  % Detect cue transitions
        cue_start = find(transitions == 1);  % Start of cue
        cue_end = find(transitions == -1) - 1;  % End of cue

        % Plot shaded regions for continuous cue periods with slight color variations
        for i = 1:length(cue_start)
            % Get the cue value and pick a variation of light orange
            cue_value = cue(cue_start(i));  
            color_idx = mod(cue_value - 1, size(color_map, 1)) + 1;  % Cycle through color variations
            fill_color = color_map(color_idx, :);  % Pick color based on cue value

            % Plot the shaded region
            fill([t(cue_start(i)), t(cue_end(i)), t(cue_end(i)), t(cue_start(i))], ...
                 [curr_ylim(1), curr_ylim(1), curr_ylim(2), curr_ylim(2)], ...
                 fill_color, 'FaceAlpha', 0.1, 'EdgeColor', 'none', 'Parent', signal_axes);
        end
    end

    % Adjust y-limits to keep the traces centered and fit the plot better
    % Centering the y-limits by setting them to [-y_offset, upper limit]
%     ylim(signal_axes, [-y_offset, (num_chan) * y_offset]);  % Adjust y-limits to better fit the traces

    % Hide y-axis ticks and labels
    set(signal_axes, 'YTick', []);  % Hide y-axis ticks
    ylabel(signal_axes, '');        % Remove y-axis label
    xlabel(signal_axes, 'Time (s)');

    % Display name of the annotated file or deafult fig title
    title(signal_axes, fig_title, 'FontSize',14, 'FontWeight', 'bold');

    % Ensure figure keeps focus after plotting
    figure(gcf);
end










