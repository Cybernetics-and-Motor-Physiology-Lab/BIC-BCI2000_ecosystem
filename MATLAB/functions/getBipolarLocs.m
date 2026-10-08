function bipolar_locs = getBipolarLocs(locs, channelLabels)
% getBipolarLocs Calculate midpoint locations for bipolar channel pairs.
%
% Inputs:
%   locs          - [nChannels x 3] electrode coordinates
%   channelLabels - channel labels in the form 'Ch%d-%d'
%
% Output:
%   bipolar_locs  - [nBipolarChannels x 3] midpoint coordinates

    bipolar_locs = nan(numel(channelLabels), 3);

    for i = 1:numel(channelLabels)

        % Extract channel indices from label, e.g. 'Ch3-4' -> [3 4]
        tokens = regexp(channelLabels{i}, 'Ch(\d+)-(\d+)', 'tokens', 'once');

        if isempty(tokens)
            warning('Could not parse bipolar channel label "%s".', channelLabels{i});
            continue;
        end

        chIdx = str2double(tokens);

        % Check that channel indices are valid
        if any(chIdx < 1) || any(chIdx > size(locs,1))
            warning('Channel indices in "%s" exceed available locations.', channelLabels{i});
            continue;
        end

        % Midpoint between the two electrode locations
        bipolar_locs(i,:) = mean(locs(chIdx,:), 1);
    end
end