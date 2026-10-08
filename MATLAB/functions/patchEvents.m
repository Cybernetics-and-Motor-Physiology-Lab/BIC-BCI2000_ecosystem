function patchEvents(ax, t, maskVector, colors)
% patchEvents Plot shaded regions corresponding to event/state masks.
%
%   patchEvents(ax, t, maskVector)
%   patchEvents(ax, t, maskVector, colors)
%
% Inputs:
%   ax         - axes handle
%   t          - time vector, same length as maskVector
%   maskVector - logical vector or numeric vector of event/state codes
%                0 = no event
%                non-zero values = event/state codes
%   colors     - optional RGB color specification:
%                [1 x 3]       -> same color for all event types
%                [N x 3]       -> one color per unique non-zero code
%
% For logical masks, true values are treated as event code 1.
%
% Example:
%   patchEvents(gca, t, artifactMask, [1 0 0])
%
%   patchEvents(gca, t, stimulusCode)
%
%   patchEvents(gca, t, stimulusCode, lines(4))

    arguments
        ax (1,1) matlab.graphics.axis.Axes
        t (:,1) double
        maskVector (:,1) {mustBeNumericOrLogical}
        colors double = []
    end

    % Convert logical mask to numeric code vector
    maskVector = double(maskVector(:));
    t = t(:);

    if isempty(maskVector)
        return
    end

    if numel(t) ~= numel(maskVector)
        error('patchEvents:SizeMismatch', ...
            't and maskVector must have the same number of elements.');
    end

    % Find unique non-zero event/state codes
    eventCodes = unique(maskVector);
    eventCodes(eventCodes == 0 | isnan(eventCodes)) = [];

    if isempty(eventCodes)
        return
    end

    nEvents = numel(eventCodes);

    % Determine colors
    if isempty(colors)
        colorMap = lines(nEvents);

    elseif isequal(size(colors), [1 3])
        % Use the same color for all event types
        colorMap = repmat(colors, nEvents, 1);

    elseif size(colors,2) == 3 && size(colors,1) == nEvents
        % One color per event code
        colorMap = colors;

    else
        error('patchEvents:InvalidColors', ...
            ['colors must be either a 1x3 RGB vector or an Nx3 matrix, ' ...
             'where N equals the number of unique non-zero event codes.']);
    end

    % Current y-axis limits
    currYLim = ylim(ax);

    % Find contiguous runs of identical mask values
    transitions = [true; diff(maskVector) ~= 0; true];

    runStart = find(transitions(1:end-1));
    runEnd   = find(transitions(2:end)) - 1;

    for i = 1:numel(runStart)

        eventValue = maskVector(runStart(i));

        % Ignore no-event periods
        if eventValue == 0 || isnan(eventValue)
            continue
        end

        % Find corresponding color
        colorIdx = find(eventCodes == eventValue, 1);
        fillColor = colorMap(colorIdx,:);

        % Start/end time
        tStart = t(runStart(i));
        tEnd   = t(runEnd(i));

        p = patch(ax, ...
            [tStart tEnd tEnd tStart], ...
            [currYLim(1) currYLim(1) currYLim(2) currYLim(2)], ...
            fillColor, ...
            'FaceAlpha', 0.15, ...
            'EdgeColor', 'none', ...
            'HandleVisibility', 'off');

        % Keep patches behind plotted signals
        uistack(p, 'bottom');
    end
end