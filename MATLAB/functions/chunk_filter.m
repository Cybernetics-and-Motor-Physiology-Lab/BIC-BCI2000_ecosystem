function filtered = chunk_filter(signal, SR, varargin)
% chunk_filter - Filter multi-channel signal with NaNs using chunked Butterworth filters
%
% Usage:
%   filtered = chunk_filter(signal, SR, 'HP', [], 'LP', [], 'notch', [], 'order', 3)
%
% Inputs:
%   signal : [samples x channels] matrix with NaNs marking gaps
%   SR     : sampling rate
%   'HP'   : high-pass cutoff frequency (Hz) or [] (default)
%   'LP'   : low-pass cutoff frequency (Hz) or [] (default)
%   'notch': center frequency to notch (Hz) and all harmonics up to Nyquist, or [] (default)
%   'order': filter order (default: 3)
%
% Output:
%   filtered : filtered signal with NaNs preserved between chunks

    % Parse inputs
    p = inputParser;
    addParameter(p, 'HP', []);
    addParameter(p, 'LP', []);
    addParameter(p, 'notch', []);
    addParameter(p, 'order', 3);
    parse(p, varargin{:});
    HP = p.Results.HP;
    LP = p.Results.LP;
    notch = p.Results.notch;
    order = p.Results.order;

    [nSamp, nChan] = size(signal);
    filtered = nan(size(signal));
    nyquist = SR / 2;

    % Design filters
    sos_all = {};  g_all = {};

    if ~isempty(HP)
        [z,p,k] = butter(order, HP/nyquist, 'high');
        [sos, g] = zp2sos(z,p,k);
        sos_all{end+1} = sos; g_all{end+1} = g;
    end

    if ~isempty(LP)
        [z,p,k] = butter(order, LP/nyquist, 'low');
        [sos, g] = zp2sos(z,p,k);
        sos_all{end+1} = sos; g_all{end+1} = g;
    end

    if ~isempty(notch)
        harmonics = notch:notch:nyquist;
        for f0 = harmonics
            bw = 2;  % Bandwidth = ±1 Hz
            Wn = [(f0 - 1), (f0 + 1)] / nyquist;
            if Wn(1) <= 0 || Wn(2) >= 1
                continue;  % skip out-of-range
            end
            [z,p,k] = butter(4, Wn, 'stop');
            [sos, g] = zp2sos(z,p,k);
            sos_all{end+1} = sos; g_all{end+1} = g;
        end
    end

    % Find continuous chunks (non-NaN regions)
    for ch = 1:nChan
        valid = ~isnan(signal(:,ch));
        edges = diff([0; valid; 0]);
        starts = find(edges == 1);
        ends = find(edges == -1) - 1;

        for i = 1:length(starts)
            idx = starts(i):ends(i);
            chunk = signal(idx, ch);

            % Skip if chunk too short
            if length(chunk) < (2 * SR)
                continue;
            end

            x = chunk;
            for f = 1:length(sos_all)
                x = filtfilt(sos_all{f}, g_all{f}, x);
            end
            filtered(idx, ch) = x;
        end
    end
end