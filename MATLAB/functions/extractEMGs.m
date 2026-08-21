function EMG_chans = extractEMGs(signal, EMG_chans, SR, fname)
% extractEMGs Compute bipolar EMG signals and annotate activity.
%
% INPUTS:
%   signal     - [samples x channels] signal matrix
%   EMG_chans  - struct with fields like:
%                   hand_idx   = [1 2]
%                   tongue_idx = [5 6]
%                   foot_idx   = [10 12]
%   SR         - sampling rate in Hz
%   fname      - filename / recording label used by pk_segment_FL
%
% OUTPUT:
%   EMG_chans  - updated struct containing fields like:
%                   EMG_chans.hand.BP_signal
%                   EMG_chans.hand.activity
%                   EMG_chans.hand.marked_idxs
%                   EMG_chans.hand.marked_times
%                   EMG_chans.hand.marked_ids
%                   EMG_chans.hand.modalities
% Dependencies: pk_segment_FL fucntion for annotation of EMG signals
    EMG_flds = fieldnames(EMG_chans);

    for f = 1:numel(EMG_flds)
        emg_idx_field = EMG_flds{f};
        % Only process fields ending in '_idx'
        if ~endsWith(emg_idx_field, '_idx')
            continue
        end
        ch_idx = EMG_chans.(emg_idx_field);
        emg_name = erase(emg_idx_field, '_idx');

        if isempty(ch_idx)
            continue
        end
        if numel(ch_idx) ~= 2
            warning('Skipping %s: expected exactly 2 channel indices.', emg_idx_field);
            continue
        end

        % Compute bipolar EMG
        EMG_bp = signal(:, ch_idx(1)) - signal(:, ch_idx(2));
        % Extract movement periods from EMG
        [marked_times, marked_ids, modalities] = ...
            pk_segment_FL(EMG_bp, [], SR, [fname emg_name]);
        % Generate activity vector
        activity = false(size(signal,1), 1);
        marked_idxs = round(marked_times .* SR);
        marked_idxs = max(1, min(marked_idxs, size(signal,1)));
        for k = 1:2:numel(marked_idxs)-1
            activity(marked_idxs(k):marked_idxs(k+1)) = true;
        end

        % Save back into structure
        EMG_chans.(emg_name).BP_signal    = EMG_bp;
        EMG_chans.(emg_name).activity     = activity;
        EMG_chans.(emg_name).marked_idxs  = marked_idxs;
        EMG_chans.(emg_name).marked_times = marked_times;
        EMG_chans.(emg_name).marked_ids   = marked_ids;
        EMG_chans.(emg_name).modalities   = modalities;
    end
end