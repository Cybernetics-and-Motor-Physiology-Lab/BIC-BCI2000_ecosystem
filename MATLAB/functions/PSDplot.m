function PSDplot(spectra,freqs,badChannels,dBAxis, transparency)
    % Validate function arguments
    arguments
        spectra {mustBeVectorOrMatrixDouble}
        freqs {mustBeVector, mustMatchLargerDimension(spectra, freqs)}
        badChannels {mustBeNumeric} = []
        dBAxis (1,1) {mustBeNumericOrLogical} = true
        transparency (1,1) = 0.2
    end

    powerLineFreq = 60; 
    % Initialize the GUI
    fig = figure('Name', 'Signal PSD', ...
            'Position', [100 100 1500 800]);
    % Create checkbox UI
    handles.checkbox = uicontrol('Style', 'checkbox', ...
        'String', 'Omit harmonics', ...
        'Units', 'normalized', ...
        'Position', [0.68 0.923 0.12 0.05], ...
        'Value', 0, ...
        'Callback', @(src, event) toggleHarmonics(src, fig, powerLineFreq));
    
    % Save original data in handles
    handles.originalSpectra = spectra;
    handles.freqs = freqs;
    handles.badChannels = badChannels;
    handles.dBAxis = dBAxis;
    handles.transparency = transparency;

    % Current spectra to plot (may be modified)
    handles.currentSpectra = spectra;

    % Plot data
    handles.ax = axes('Parent', fig);
    handles.plotHandles = plotPSD(handles);

    % Store all data in figure
    set(fig, 'UserData', handles);
end
%% Helper Functions

function plotHandles = plotPSD(handles)
    cla(handles.ax); % Clear previous plot
    hold(handles.ax, 'on');

    spectra = handles.currentSpectra;
    freqs = handles.freqs;
    badChannels = handles.badChannels;
    dBAxis = handles.dBAxis;
    transparency = handles.transparency;
    
    chnls = 1:size(spectra,1);
    chnls(badChannels) = [];     % Get rid of bad channels
    
    plotHandles = gobjects(length(chnls) + 1, 1);

   if dBAxis
        for ch = 1:length(chnls)
            plotHandles(ch) = plot(handles.ax, freqs, pow2db(spectra(chnls(ch),:)), ...
                'color', [0 0.4470 0.7410 transparency], 'LineWidth', 1.5);
            plotHandles(ch).DataTipTemplate.DataTipRows(end) = ...
                dataTipTextRow("Channel", repelem(chnls(ch), numel(freqs)));
        end
        plotHandles(end) = plot(handles.ax, freqs, pow2db(nanmean(spectra(chnls,:),1)), ...
            'color', [0.8500 0.3250 0.0980], 'DisplayName', 'Average', 'LineWidth', 2.5);
        ylabel(handles.ax, 'Power (dB/Hz)');
    else
        for ch = 1:length(chnls)
            plotHandles(ch) = semilogy(handles.ax, freqs, spectra(chnls(ch),:), ...
                'color', [0 0.4470 0.7410 transparency], 'LineWidth', 1.5);
            plotHandles(ch).DataTipTemplate.DataTipRows(end) = ...
                dataTipTextRow("Channel", repelem(chnls(ch), numel(freqs)));
        end
        plotHandles(end) = semilogy(handles.ax, freqs, nanmean(spectra(chnls,:),1), ...
            'color', [0.8500 0.3250 0.0980], 'DisplayName', 'Average', 'LineWidth', 2.5);
        ylabel(handles.ax, 'Power (\muV^2/Hz)');
   end

    title(handles.ax, 'Noise floor');
    xlim(handles.ax, [freqs(1), freqs(end)]);
    xlabel(handles.ax, 'Frequency (Hz)');
    legend(legendUnq(handles.ax));
    grid(handles.ax, 'on');
    fontsize(20,'points');
end

function toggleHarmonics(src, fig, plFreq)
    handles = get(fig, 'UserData');
    stopWidth = 2;
    if get(src, 'Value') == 1  % Checkbox checked
        spectra = handles.originalSpectra;
        freqs = handles.freqs;

        mask = false(size(freqs));
        maxF = max(freqs);
        
        for k = 1:floor(maxF/plFreq)
            center = k * plFreq;
            mask(freqs >= center-stopWidth & freqs <= center+stopWidth) = true;
        end
        spectra(:, mask) = NaN;
        handles.currentSpectra = spectra;
    else
        handles.currentSpectra = handles.originalSpectra;
    end

    % Redraw plot
    handles.plotHandles = plotPSD(handles);

    % Store updated data
    set(fig, 'UserData', handles);
end

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
