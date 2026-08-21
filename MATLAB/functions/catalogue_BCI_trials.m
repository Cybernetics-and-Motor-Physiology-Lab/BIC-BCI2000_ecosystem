function trialInfo = catalogue_BCI_trials(states, fs, postTrialWindow, verbose)
% catalogue_BCI_trials Catalogue BCI2000 CursorTask trials from state vectors.
% Extracts trial starts and ends based on the Feedback and TargetCode
%
% Syntax
%   trialInfo = catalogue_BCI_trials(states)
%   trialInfo = catalogue_BCI_trials(states, fs)
%   trialInfo = catalogue_BCI_trials(states, fs, postTrialWindow)
%   trialInfo = catalogue_BCI_trials(states, fs, postTrialWindow, verbose)
%
% Inputs
%   states           - struct with fields:
%                        states.TargetCode
%                        states.ResultCode
%                        states.Feedback
%   fs               - sampling rate (Hz), required if postTrialWindow > 0
%   postTrialWindow  - seconds after trial end to inspect ResultCode (default = 0)
%   verbose          - logical flag to print summary (default = false)
%
% Output
%   trialInfo - struct with fields:
%       .target
%       .performance
%       .startIdx
%       .endIdx
%       .nTrials
%       .nHits
%       .nMisses
%       .nInvalid
%       .summary

    % Defaults
    if nargin < 2 || isempty(fs), fs = []; end
    if nargin < 3 || isempty(postTrialWindow), postTrialWindow = 0; end
    if nargin < 4 || isempty(verbose), verbose = false; end

    % Validation
    if ~isstruct(states)
        error('states must be a struct.');
    end
    if ~isfield(states, 'TargetCode') || ~isfield(states, 'ResultCode') || ~isfield(states, 'Feedback')
        error('states must contain TargetCode, ResultCode and Feedback.');
    end
    if postTrialWindow > 0 && isempty(fs)
        error('fs must be provided when postTrialWindow > 0.');
    end

    target = double(states.TargetCode);
    result = double(states.ResultCode);
    feedback = double(states.Feedback);

    nSamples = numel(target);
    BCIcontrol = zeros(nSamples,1);

    % --- Trial detection based on Target Code---
    % Start: 0 -> nonzero
    targetStart = find(target(1:end-1) == 0 & target(2:end) > 0) + 1;
   
    if target(1) > 0
        targetStart = [1; targetStart];
    end

    % End: nonzero -> 0
    targetEnd = find(target(1:end-1) > 0 & target(2:end) == 0);
    if target(end) > 0
        targetEnd = [targetEnd; nSamples];
    end

    nTrials = numel(targetStart);
    trialTarget = target(targetStart);
    trialPerf = nan(nTrials,1);
    trialStart = nan(nTrials,1);
    trialEnd = nan(nTrials,1);

    postSamples = round(postTrialWindow * fs);

    % Iterate over trials and extract feedback times and results
    for k = 1:nTrials
        idx1 = targetStart(k);
        idx2_target = targetEnd(k);

        idx2_result = idx2_target;
        if postSamples > 0
            idx2_result = min(idx2_target + postSamples, nSamples);
        end

        % Window for feedback/control period only
        idx_feedback = idx1:idx2_target;

        % Window for result inspection, optionally extended
        idx_result = idx1:idx2_result;

        trialFeedback = feedback(idx_feedback);
        trialResult = result(idx_result);
        thisTarget = trialTarget(k);

        % Feedback onset: 0 -> 1 within this target period
        fbStartLocal = find(trialFeedback(1:end-1) == 0 & trialFeedback(2:end) > 0, ...
            1, 'first') + 1;

        % Feedback offset: 1 -> 0 within this target period
        fbEndLocal = find(trialFeedback(1:end-1) > 0 & trialFeedback(2:end) == 0, ...
            1, 'last');

        % Handle case where feedback is already on at idx1
        if isempty(fbStartLocal) && trialFeedback(1) > 0
            fbStartLocal = 1;
        end

        % Handle case where feedback stays on until target end
        if isempty(fbEndLocal) && any(trialFeedback > 0)
            fbEndLocal = numel(trialFeedback);
        end

        % Convert local feedback indices to global indices
        if ~isempty(fbStartLocal) && ~isempty(fbEndLocal) && fbEndLocal >= fbStartLocal
            trialStart(k) = idx_feedback(fbStartLocal);
            trialEnd(k)   = idx_feedback(fbEndLocal);

            BCIcontrol(trialStart(k):trialEnd(k)) = thisTarget;
        else
            trialStart(k) = NaN;
            trialEnd(k)   = NaN;
        end

        % Extract results
        if any(trialResult == thisTarget)
            trialPerf(k) = 1;
        elseif any(trialResult > 0)
            trialPerf(k) = -1;
        else
            trialPerf(k) = 0;
        end
    end

    nHits    = sum(trialPerf == 1);
    nMisses  = sum(trialPerf == -1);
    nInvalid = sum(trialPerf == 0);

    summaryStr = ['Total targets shown: ' num2str(nTrials) ...
        ', ' num2str(nHits) ' correct' ...
        ', ' num2str(nInvalid) ' invalid/timeouts' ...
        ', ' num2str(nMisses) ' incorrect'];

    % Verbose output
    if verbose
        disp(summaryStr);
    end

    trialInfo = struct( ...
        'targetCode', trialTarget, ...
        'BCIcontrol', BCIcontrol, ...
        'performance', trialPerf, ...
        'startIdx', trialStart, ...
        'endIdx', trialEnd, ...
        'nTrials', nTrials, ...
        'nHits', nHits, ...
        'nMisses', nMisses, ...
        'nInvalid', nInvalid, ...
        'summary', summaryStr);
end