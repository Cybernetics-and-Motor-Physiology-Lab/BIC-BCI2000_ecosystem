function spectStats = trialSpectralStats(trialSpectralFeatures, LFB_range, HFB_range, BF_correct)
% trialSpectralStats Compute trial-wise spectral statistics.
%
% Inputs
%   trialSpectralFeatures - struct containing trial spectral features.
%   LFB_range             - [fLow fHigh] low-frequency band in Hz.
%   HFB_range             - [fLow fHigh] high-frequency band in Hz.
%   BF_correct            - optional logical, Bonferroni correction. Default true.
%
% Output
%   spectStats            - struct containing spectral statistics.

    if nargin < 4 || isempty(BF_correct)
        BF_correct = true;
    end

    %% Validate input structure
    requiredFields = {'trialNormPSD','freqBins','validTrial','trialStimCode', ...
        'trialStartIdx','trialStopIdx','trialStartTime','trialStopTime'};

    if ~isstruct(trialSpectralFeatures) || isempty(trialSpectralFeatures)
        error('trialSpectralFeatures must be a nonempty struct.');
    end

    for i = 1:numel(requiredFields)
        if ~isfield(trialSpectralFeatures, requiredFields{i})
            error('trialSpectralFeatures is missing required field "%s".', requiredFields{i});
        end
    end

    validateattributes(LFB_range, {'numeric'}, {'vector','finite','increasing'});
    validateattributes(HFB_range, {'numeric'}, {'vector','finite','increasing'});
    validateattributes(BF_correct, {'logical','numeric'}, {'scalar'});

    %% Remove invalid trials
    validMask = logical(trialSpectralFeatures.validTrial(:));
    bad_trials = find(~validMask);

    for i = 1:numel(bad_trials)
        fprintf('Ignoring trial %d\n', bad_trials(i));
    end

    tsf = removeInvalidTrials(trialSpectralFeatures, validMask);

    nps = tsf.trialNormPSD;              % [freq x channels x trials]
    f = tsf.freqBins(:).';
    tr_sc = tsf.trialStimCode(:).';      % [1 x trials]

    num_chans = size(nps,2);
    nFreq = numel(f);

    if size(nps,1) ~= nFreq
        error('Dimension mismatch: trialNormPSD second dimension must match freqBins.');
    end

    if size(nps,3) ~= numel(tr_sc)
        error('Dimension mismatch: trialNormPSD third dimension must match trialStimCode.');
    end

    %% Find band indices
    LFB_idxs = ismember(round(f,6), round(LFB_range(:),6));
    HFB_idxs = ismember(round(f,6), round(HFB_range(:),6));

    if ~any(LFB_idxs)
        error('No frequency bins found inside LFB_range [%g %g].', LFB_range(1), LFB_range(2));
    end

    if ~any(HFB_idxs)
        error('No frequency bins found inside HFB_range [%g %g].', HFB_range(1), HFB_range(2));
    end

    %% Trial groups
    taskMask = tr_sc == 1;
    restMask = tr_sc == 0;

    if ~any(taskMask)
        warning('No task trials found where trialStimCode == 1.');
    end

    if ~any(restMask)
        warning('No rest trials found where trialStimCode == 0.');
    end

    %% Legacy PSD/r2 statistics
    HFB_trials = squeeze(mean(nps(HFB_idxs,:, :), 1, 'omitnan')); % [ch x trials]
    LFB_trials = squeeze(mean(nps(LFB_idxs,:, :), 1, 'omitnan')); % [ch x trials]

    if num_chans == 1
        HFB_trials = reshape(HFB_trials, 1, []);
        LFB_trials = reshape(LFB_trials, 1, []);
    end

    r_HFB = nan(1, num_chans);
    r_LFB = nan(1, num_chans);
    p_HFB = nan(1, num_chans);
    p_LFB = nan(1, num_chans);
    rmap = nan(num_chans, nFreq);

    for chan = 1:num_chans
        r_HFB(chan) = rsa(HFB_trials(chan, taskMask), HFB_trials(chan, restMask));
        r_LFB(chan) = rsa(LFB_trials(chan, taskMask), LFB_trials(chan, restMask));

        [~, p_HFB(chan)] = ttest2(HFB_trials(chan, taskMask), HFB_trials(chan, restMask));
        [~, p_LFB(chan)] = ttest2(LFB_trials(chan, taskMask), LFB_trials(chan, restMask));

        if BF_correct
            p_HFB(chan) = min(p_HFB(chan) * num_chans, 1);
            p_LFB(chan) = min(p_LFB(chan) * num_chans, 1);
        end

        for fi = 1:nFreq
            rmap(chan, fi) = rsa(squeeze(nps(fi, chan, taskMask)), ...
                squeeze(nps(fi, chan, restMask)));
        end
    end

    %% Integrated power in selected bands from normalized PSD
    HFB_power = squeeze(trapz(f(HFB_idxs), nps(HFB_idxs,:,:), 1));
    LFB_power = squeeze(trapz(f(LFB_idxs), nps(LFB_idxs,:,:), 1));

    if num_chans == 1
        HFB_power = reshape(HFB_power, 1, []);
        LFB_power = reshape(LFB_power, 1, []);
    end

    bandStats = struct();
    bandStats.HFB.power = HFB_power;
    bandStats.LFB.power = LFB_power;

    bandStats.HFB.task = HFB_power(:, taskMask);
    bandStats.HFB.rest = HFB_power(:, restMask);
    bandStats.LFB.task = LFB_power(:, taskMask);
    bandStats.LFB.rest = LFB_power(:, restMask);

    bandStats.HFB.r = nan(1,num_chans);
    bandStats.HFB.p = nan(1,num_chans);
    bandStats.LFB.r = nan(1,num_chans);
    bandStats.LFB.p = nan(1,num_chans);

    for chan = 1:num_chans
        bandStats.HFB.r(chan) = rsa(bandStats.HFB.task(chan,:), bandStats.HFB.rest(chan,:));
        bandStats.LFB.r(chan) = rsa(bandStats.LFB.task(chan,:), bandStats.LFB.rest(chan,:));

        [~, bandStats.HFB.p(chan)] = ttest2(bandStats.HFB.task(chan,:), bandStats.HFB.rest(chan,:));
        [~, bandStats.LFB.p(chan)] = ttest2(bandStats.LFB.task(chan,:), bandStats.LFB.rest(chan,:));

        if BF_correct
            bandStats.HFB.p(chan) = min(bandStats.HFB.p(chan) * num_chans, 1);
            bandStats.LFB.p(chan) = min(bandStats.LFB.p(chan) * num_chans, 1);
        end
    end

    %% Statistics for precomputed bandPower fields, if available
    precomputedBandStats = struct();

    if isfield(tsf, 'bandPower') && isstruct(tsf.bandPower)
        bpNames = fieldnames(tsf.bandPower);

        for b = 1:numel(bpNames)
            bandName = bpNames{b};
            bp = tsf.bandPower.(bandName); % expected [channels x trials]

            if size(bp,2) ~= numel(tr_sc)
                warning('Skipping bandPower.%s due to trial dimension mismatch.', bandName);
                continue
            end

            precomputedBandStats.(bandName).task = bp(:, taskMask);
            precomputedBandStats.(bandName).rest = bp(:, restMask);
            precomputedBandStats.(bandName).r = nan(1,num_chans);
            precomputedBandStats.(bandName).p = nan(1,num_chans);

            for chan = 1:num_chans
                precomputedBandStats.(bandName).r(chan) = rsa( ...
                    bp(chan, taskMask), bp(chan, restMask));

                [~, precomputedBandStats.(bandName).p(chan)] = ttest2( ...
                    bp(chan, taskMask), bp(chan, restMask));

                if BF_correct
                    precomputedBandStats.(bandName).p(chan) = ...
                        min(precomputedBandStats.(bandName).p(chan) * num_chans, 1);
                end
            end
        end
    end

    %% Coherence statistics
    coherenceStats = struct();

    if isfield(tsf, 'trialCoherence') && ~isempty(tsf.trialCoherence)
        cohTask = tsf.trialCoherence(taskMask);
        cohRest = tsf.trialCoherence(restMask);

        cohTask = cohTask(~cellfun(@isempty, cohTask));
        cohRest = cohRest(~cellfun(@isempty, cohRest));

        if ~isempty(cohTask)
            coherenceStats.meanTaskCoherence = mean(cat(4, cohTask{:}), 4, 'omitnan');
            coherenceStats.rmsTaskCoherence = sqrt(mean(coherenceStats.meanTaskCoherence.^2, 3, 'omitnan'));
        else
            coherenceStats.meanTaskCoherence = [];
            coherenceStats.rmsTaskCoherence = [];
        end

        if ~isempty(cohRest)
            coherenceStats.meanRestCoherence = mean(cat(4, cohRest{:}), 4, 'omitnan');
            coherenceStats.rmsRestCoherence = sqrt(mean(coherenceStats.meanRestCoherence.^2, 3, 'omitnan'));
        else
            coherenceStats.meanRestCoherence = [];
            coherenceStats.rmsRestCoherence = [];
        end

        if ~isempty(cohTask) && ~isempty(cohRest)
            coherenceStats.deltaMeanCoherence = ...
                coherenceStats.meanTaskCoherence - coherenceStats.meanRestCoherence;

            coherenceStats.deltaRMSCoherence = ...
                coherenceStats.rmsTaskCoherence - coherenceStats.rmsRestCoherence;
        else
            coherenceStats.deltaMeanCoherence = [];
            coherenceStats.deltaRMSCoherence = [];
        end
    else
        coherenceStats.meanTaskCoherence = [];
        coherenceStats.meanRestCoherence = [];
        coherenceStats.rmsTaskCoherence = [];
        coherenceStats.rmsRestCoherence = [];
        coherenceStats.deltaMeanCoherence = [];
        coherenceStats.deltaRMSCoherence = [];
    end

    %% Wrap results
    spectStats = struct();

    spectStats.rvals = struct();
    spectStats.rvals.HFB_range = HFB_range;
    spectStats.rvals.LFB_range = LFB_range;
    spectStats.rvals.HFB_idxs = HFB_idxs;
    spectStats.rvals.LFB_idxs = LFB_idxs;
    spectStats.rvals.HFB_trials = HFB_trials;
    spectStats.rvals.LFB_trials = LFB_trials;
    spectStats.rvals.r_HFB = r_HFB;
    spectStats.rvals.r_LFB = r_LFB;
    spectStats.rvals.p_HFB = p_HFB;
    spectStats.rvals.p_LFB = p_LFB;
    spectStats.rvals.rmap = rmap;

    spectStats.bandStats = bandStats;
    spectStats.precomputedBandStats = precomputedBandStats;
    spectStats.coherenceStats = coherenceStats;

    spectStats.trialInfo = struct();
    spectStats.trialInfo.validMaskOriginal = validMask;
    spectStats.trialInfo.removedTrials = bad_trials;
    spectStats.trialInfo.trialStimCode = tr_sc;
    spectStats.trialInfo.taskMask = taskMask;
    spectStats.trialInfo.restMask = restMask;

    spectStats.info = struct();
    spectStats.info.BF_correct = logical(BF_correct);
    spectStats.info.freqBins = f;
end


function tsf = removeInvalidTrials(tsf, validMask)

    trialFields3D = {'trialPSD','trialNormPSD'};
    for i = 1:numel(trialFields3D)
        fld = trialFields3D{i};
        if isfield(tsf, fld) && ~isempty(tsf.(fld))
            tsf.(fld) = tsf.(fld)(:,:,validMask);
        end
    end

    trialFieldsVec = {'validTrial','trialStimCode','trialStartIdx','trialStopIdx', ...
        'trialStartTime','trialStopTime'};
    for i = 1:numel(trialFieldsVec)
        fld = trialFieldsVec{i};
        if isfield(tsf, fld) && ~isempty(tsf.(fld))
            tsf.(fld) = tsf.(fld)(validMask);
        end
    end

    if isfield(tsf, 'bandPower') && isstruct(tsf.bandPower)
        bpNames = fieldnames(tsf.bandPower);
        for b = 1:numel(bpNames)
            fld = bpNames{b};
            if size(tsf.bandPower.(fld),2) == numel(validMask)
                tsf.bandPower.(fld) = tsf.bandPower.(fld)(:,validMask);
            end
        end
    end

    if isfield(tsf, 'trialCPSD') && ~isempty(tsf.trialCPSD)
        tsf.trialCPSD = tsf.trialCPSD(validMask);
    end

    if isfield(tsf, 'trialCoherence') && ~isempty(tsf.trialCoherence)
        tsf.trialCoherence = tsf.trialCoherence(validMask);
    end
end

%% R2 activation helper 

function outrsa=rsa(d1,d2)
% function outrsa=rsa(d1,d2)
% this function calculates the signed r-squared cross-correlation for
% vectors d1 and d2.  it is signed to reflect d1>d2
% kjm 2007

d1=reshape(d1,1,[]);
d2=reshape(d2,1,[]);
outrsa=((mean(d1)-mean(d2))^3)/abs(mean(d1)-mean(d2))...
    /var([d1 d2])...
    *(length(d1)*length(d2))/length([d1 d2])^2;
end