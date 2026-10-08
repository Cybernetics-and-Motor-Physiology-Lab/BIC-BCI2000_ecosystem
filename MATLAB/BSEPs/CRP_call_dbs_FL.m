function [crp_projs, crp_parms, bad_trials, baselines] = CRP_call_dbs_FL( ...
    data, t, opts)
    % CRP_call_dbs_FL
    %
    % Flexible wrapper for CRP_method with exposed parameters.
    %
    % USAGE:
    %   [crp_projs, crp_parms] = CRP_call_dbs_FL(data, t)
    %   [crp_projs, crp_parms] = CRP_call_dbs_FL(data, t, opts)
    %
    % opts fields (all optional, defaults shown):
    %   opts.srate            = 4800;
    %   opts.stim_time        = 0;        % seconds
    %   opts.t1               = 0.005;    % post-stim start (s)
    %   opts.t2               = 1.0;      % post-stim end (s)
    %   opts.t_baseline1      = 0.25;     % baseline start before stim (s)
    %   opts.t_baseline2      = 0.005;    % baseline end before stim (s)
    %   opts.p_thresh         = 1e-6;
    %   opts.artifact_removal = true;
    %   opts.art_interval     = 'tR';     % 'tR' or 'full'
    %
    % kjm 7/2022
    % Refactored FL 01/2026
    
    %% -------------------- Defaults --------------------
    if nargin < 3 || isempty(opts)
        opts = struct();
    end
    
    opts = setDefault(opts, 'srate',            4800);
    opts = setDefault(opts, 'stim_time',        0);
    opts = setDefault(opts, 't1',               0.005);
    opts = setDefault(opts, 't2',               1.0);
    opts = setDefault(opts, 't_baseline1',      0.25);
    opts = setDefault(opts, 't_baseline2',      0.005);
    opts = setDefault(opts, 'p_thresh',          1e-6);
    opts = setDefault(opts, 'artifact_removal',  true);
    opts = setDefault(opts, 'art_interval',     'tR');
    
    %% -------------------- Derived parameters --------------------
    srate = opts.srate;
    
    stim_idx = find(t == opts.stim_time, 1, 'first');
    if isempty(stim_idx)
        error('CRP_call_dbs_FL:StimNotFound', ...
            'No sample found at stim_time = %.6f s', opts.stim_time);
    end
    
    art_rem.do       = opts.artifact_removal; 
    art_rem.interval = opts.art_interval;
    
    %% -------------------- Time window selection --------------------
    if t(1) >= opts.t1
        error('CRP_call_dbs_FL:BadTimeLimits', ...
            't(1) >= t1 (%.4f s)', opts.t1);
    elseif t(end) < opts.t2
        error('CRP_call_dbs_FL:BadTimeLimits', ...
            't(end) < t2 (%.4f s)', opts.t2);
    end
    
    tpts = find(t > opts.t1 & t <= opts.t2);
    t_win = t(tpts);
    
    %% -------------------- Baseline subtraction --------------------
    num_trials = size(data, 2);
    baselines  = zeros(num_trials, 1);
    
    b1 = floor(stim_idx - opts.t_baseline1 * srate);
    b2 = floor(stim_idx - opts.t_baseline2 * srate);
    
    if b1 < 1 || b2 < 1 || b1 >= b2
        error('CRP_call_dbs_FL:BadBaselineWindow', ...
            'Invalid baseline window relative to stim index.');
    end
    
    for k = 1:num_trials
        baselines(k) = mean(data(b1:b2, k), 'omitnan');
        data(:, k)   = data(:, k) - baselines(k);
    end
    
    V = data(tpts, :);
    
    %% -------------------- Artifact rejection --------------------
    bad_trials = NaN;
    
    if art_rem.do
        [~, crp_tmp] = CRP_method(V, t_win);
        [m, p] = trial_tests(crp_tmp, art_rem, size(V, 2));
        % conditions to be potentially artifactual 
        % 1 -- p<threshold 
        % 2 -- represents a decrease in projection magnitude
        bad_idx = find(p < opts.p_thresh & m < mean(m));
    
        if ~isempty(bad_idx)
            for k = bad_idx(:)'
                disp(['Trial ' num2str(k) ...
                    ' appears artifactual, p=' num2str(p(k))])
            end
            disp(['Automatically rejecting trials: ' num2str(bad_idx)])
            V(:, bad_idx) = [];
        end
    
        bad_trials = bad_idx;
    end
    
    %% -------------------- Final CRP call --------------------
    [crp_parms, crp_projs] = CRP_method(V, t_win);

end

function [m, p] = trial_tests(crp_projs, art_rem, num_trials)
    % [m,p]=trial_tests(crp_projs,art_rem,num_trials)
    % Identifies outliers in individual trials
    % kjm, 7/2022
    
    %% Select test interval (unchanged behavior)
    if strcmp(art_rem.interval, 'tR')
        S_test = crp_projs.S_all(:, crp_projs.tR_index);
    elseif strcmp(art_rem.interval, 'full')
        S_test = crp_projs.S_all(:, end);
    else
        error('trial_tests:BadInterval', ...
            'art_rem.interval must be ''tR'' or ''full''.');
    end
    
    %% Preallocate outputs (same results, faster/cleaner)
    m = zeros(1, num_trials);
    p = zeros(1, num_trials);
    
    %% Account for double counting (logic preserved)
    stat_indices = 1:2:(num_trials^2 - num_trials);
    
    if rem(num_trials, 2) == 1  % odd number of trials
        b = zeros(size(stat_indices)); % offset mask
        halfBlock = (num_trials - 1) / 2;
    
        for k = 1:num_trials
            if mod(k, 2) == 0
                idx1 = (k-1)*halfBlock + 1;
                idx2 = k*halfBlock;
                b(idx1:idx2) = 1;
            end
        end
    
        stat_indices = stat_indices + b;
    end
    
    excl_indices = setdiff(1:numel(S_test), stat_indices);
    % End of part for double counting
    
    %% Main loop (same indices, same stats)
    for q = 1:num_trials
    
        if q > 1
            t_indices1 = (q-1) : (num_trials-1) : (q*(num_trials-1));
        else
            t_indices1 = [];
        end
    
        t_indices2 = (q*num_trials) : (num_trials-1) : numel(S_test);
        t_indices3 = (q-1)*(num_trials-1) + (1:(num_trials-1));
    
        g_in  = unique([t_indices1, t_indices2, t_indices3]);
        g_out = setdiff(1:numel(S_test), [g_in, excl_indices]);
    
        [~, p(q)] = ttest2(S_test(g_in), S_test(g_out));
        m(q)      = mean(S_test(t_indices3));
    
    end
end

function s = setDefault(s, field, value)
    if ~isfield(s, field) || isempty(s.(field))
        s.(field) = value;
    end
end



