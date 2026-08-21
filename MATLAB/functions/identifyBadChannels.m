function [badChannels, badChannelInfo] = identifyBadChannels(data, construct, FS, k, knownBadChnls, impedanceValues, graphicOut)
% identifyBadChannels Detect bad channels in multichannel electrophysiology data.
%
% Syntax
%   badChannels = identifyBadChannels(data, construct, FS, k)
%   [badChannels, badChannelInfo] = identifyBadChannels(data, construct, FS, k, knownBadChnls, impedanceValues, graphicOut)
%
% Inputs
%   data             - [samples x channels] numeric signal matrix.
%   construct        - Electrode construct:
%                        'surface'
%                        'surface_wide'
%                        'depth'
%                        'combined'
%   FS               - Sampling rate in Hz (mandatory).
%   k                - Robust threshold multiplier (number of "robust sigmas"
%                      from the median). Used in both time and spectral domain.
%                      Default = 3.
%   knownBadChnls    - Optional vector of known bad channel indices.
%                      These are directly marked as bad and excluded from
%                      all subsequent statistical estimation.
%                      Default = [].
%   impedanceValues  - Optional Nx2 cell array:
%                      Column 1: channel label or index (e.g. 'Ch12' or 12)
%                      Column 2: impedance value
%                      Channels with impedance > 10000 Ohm or Inf are marked bad.
%                      Default = [].
%   graphicOut       - Optional logical flag to display diagnostic plots.
%                      Default = false.
%
% Outputs
%   badChannels      - Numeric vector of detected bad channel indices.
%   badChannelInfo   - Struct containing:
%                        .badMask        - Logical vector marking bad channels
%                        .metrics        - Structure with computed features:
%                                             .channelSTD
%                                             .spectralSTD
%                                             .eegPower
%                                             .noiseFloor
%                        .thresholds     - Structure with thresholds:
%                                             .globalTimeThreshold
%                                             .cleanTimeThreshold
%                                             .leadThresholds
%                                             .k
%                        .stats          - Structure with diagnostic statistics:
%                                             .zEEG
%                                             .zNoise
%                                             .outEEG
%                                             .outNoise
%
% Method overview
%   The algorithm detects bad channels using a multi-stage robust pipeline:
%   1. Prior exclusion
%       - Channels listed in knownBadChnls are marked as bad.
%       - Channels with high impedance (>10 kOhm or Inf) are marked as bad.
%   2. Time-domain detection
%       - Compute channel-wise standard deviation.
%       - Estimate global robust statistics (median + MAD).
%       - Identify candidate channels exceeding k * robust sigma.
%   3. Spectral-domain detection
%       - Compute Welch PSD.
%       - Extract:
%           • EEG-band power (1 Hz → 150Hz)
%           • High-frequency noise floor (0.75*Nyquist → Nyquist)
%       - Estimate spectral variability across channels.
%   4. Lead-wise refinement
%       - Compute robust thresholds within each lead.
%       - Adjust thresholds if lead variability is elevated.
%   5. Candidate aggregation
%       - Combine time-domain, spectral, and lead-based candidates.
%   6. Robust spectral confirmation
%       - Compute robust z-scores:
%             z = |x - median| / (1.4826 * MAD)
%       - Detect outliers using:
%             • robust z-score threshold (z > k)
%             • MATLAB isoutlier(...,'median')
%       - Channels exceeding thresholds in EEG power or noise floor
%         are confirmed as bad.
%   7. Final decision
%       - Combine:
%             known bad + impedance bad + confirmed candidates
%       - Return final bad channel indices.
%
% Notes
%   - The method avoids parametric statistical tests (e.g., t-tests) and
%     instead uses robust statistics suited for small channel counts.
%   - All thresholds are based on median and MAD, making the method robust
%     to outliers.
%   - Spectral estimation uses adaptive windowing when signals are short.
%
% Author
%   Frederik Lampert, Mayo Clinic, 2025

    %% -------------------- Defaults and validation --------------------
    if nargin < 4 || isempty(k)
        k = 3;
    end

    if nargin < 5 || isempty(knownBadChnls)
        knownBadChnls = [];
    end

    if nargin < 6
        impedanceValues = [];
    end

    if nargin < 7 || isempty(graphicOut)
        graphicOut = false;
    end

    validateattributes(data, {'numeric'}, {'2d','nonempty','real'}, mfilename, 'data', 1);
    validateattributes(FS, {'numeric'}, {'scalar','positive','finite'}, mfilename, 'FS', 3);
    validateattributes(k, {'numeric'}, {'scalar','positive','finite'}, mfilename, 'k', 4);
    validateattributes(graphicOut, {'logical','numeric'}, {'scalar'}, mfilename, 'graphicOut', 7);

    construct = string(construct);
    validConstructs = ["surface","surface_wide","depth","combined"];
    if ~ismember(construct, validConstructs)
        error('Invalid construct "%s". Must be one of: %s.', construct, strjoin(validConstructs, ', '));
    end

    nSamp = size(data,1);
    nCh = size(data,2);

    knownBadChnls = knownBadChnls(:).';
    if any(knownBadChnls < 1 | knownBadChnls > nCh)
        error('knownBadChnls contains indices outside the valid channel range 1:%d.', nCh);
    end

    %% -------------------- Define leads --------------------
    switch construct
        case "surface"
            leads = {1:12; 13:22; 23:32};

        case "surface_wide"
            leads = {1:16; 17:32};

        case "depth"
            leads = {1:8; 9:16; 17:24; 25:32};

        case "combined"
            leads = {1:8; 9:16; 17:24; 25:32};
    end

    % Clip lead definitions to existing channel count
    for l = 1:numel(leads)
        leads{l} = leads{l}(leads{l} >= 1 & leads{l} <= nCh);
    end

    %% -------------------- Initialize masks --------------------
    knownBadMask = false(1,nCh);
    knownBadMask(knownBadChnls) = true;
    impedanceBadMask = false(1,nCh);
    timeCandidateMask = false(1,nCh);
    spectralCandidateMask = false(1,nCh);
    leadCandidateMask = false(1,nCh);
    confirmedSpectralMask = false(1,nCh);
    badMask = false(1,nCh);

    %% -------------------- Impedance-based bad channels --------------------
    if ~isempty(impedanceValues)
        if ~iscell(impedanceValues) || size(impedanceValues,2) < 2
            warning('impedanceValues must be an Nx2 cell array. Ignoring impedance values.');
        else
            for r = 1:size(impedanceValues,1)
                chIdx = parseChannelIndex(impedanceValues{r,1});
                impVal = impedanceValues{r,2};

                if iscell(impVal)
                    impVal = impVal{1};
                end
                if ischar(impVal) || isstring(impVal)
                    impVal = str2double(string(impVal));
                end

                if ~isempty(chIdx) && chIdx >= 1 && chIdx <= nCh
                    if isnumeric(impVal) && isscalar(impVal)
                        if isinf(impVal) || impVal > 10000
                            impedanceBadMask(chIdx) = true;
                        end
                    end
                end
            end
        end
    end

    badMask = badMask | knownBadMask | impedanceBadMask;

    %% -------------------- Analysis mask --------------------
    analysisMask = ~badMask;

    if sum(analysisMask) < 2
        warning('Fewer than two channels available for bad-channel estimation after known/impedance exclusions.');
        badChannels = find(badMask);
        badChannelInfo = packInfo();
        return
    end

    %% -------------------- Time-domain robust channel STD --------------------
    chanSTD = std(data, 0, 1, 'omitnan');

    globalTimeMedian = median(chanSTD(analysisMask), 'omitnan');
    globalTimeSigma = robustSigmaMAD(chanSTD(analysisMask));
    globalTimeThresh = globalTimeMedian + k * globalTimeSigma;

    timeCandidateMask = chanSTD > globalTimeThresh;
    timeCandidateMask = timeCandidateMask & analysisMask;

    %% -------------------- Welch PSD --------------------
    winLength = min(round(10 * FS), nSamp);

    if winLength < round(1 * FS)
        warning('Signal is shorter than 1 second. Spectral bad-channel estimation may be unreliable.');
        freqVec = [];
        pxx = [];
        logPxx = [];
        eegPower = nan(1,nCh);
        noiseFloor = nan(1,nCh);
        spectralSTD = nan(1,nCh);
    else
        winLength = max(winLength, min(nSamp, round(1 * FS)));
        wnd = hann(winLength, 'periodic');
        noverlap = floor(0.5 * winLength);

        % Approximate 1 Hz bins where possible
        nfft = max(2^nextpow2(winLength), 2^nextpow2(round(FS)));

        [pxx, freqVec] = pwelch(data, wnd, noverlap, nfft, FS);
        % pxx: [freq x channels]

        logPxx = log10(pxx + eps);

        nyq = FS / 2;
        % Specify the EEG bandas 1-150Hz
        eegBand = freqVec >= 1 & freqVec <= min(150, nyq);
        % Specify the noise-floor frequency range
        noiseLow = max(250, 0.75*nyq);
        noiseBand = freqVec >= noiseLow & freqVec <= nyq;

        if ~any(noiseBand)
            warning('Noise band is empty for FS = %.1f Hz. Noise-floor metric set to NaN.', FS);
            noiseFloor = nan(1, nCh);
        else
            noiseFloor = mean(logPxx(noiseBand,:), 1, 'omitnan');
        end

        eegPower = mean(logPxx(eegBand,:), 1, 'omitnan');

        spectralSTD = std(logPxx, 0, 1, 'omitnan');

        globalSpecMedian = median(spectralSTD(analysisMask), 'omitnan');
        globalSpecSigma = robustSigmaMAD(spectralSTD(analysisMask));
        globalSpecThresh = globalSpecMedian + k * globalSpecSigma;

        spectralCandidateMask = spectralSTD > globalSpecThresh;
        spectralCandidateMask = spectralCandidateMask & analysisMask;
    end

    %% -------------------- Recalculate global stats after first candidates --------------------
    firstCandidateMask = timeCandidateMask | spectralCandidateMask;
    cleanedAnalysisMask = analysisMask & ~firstCandidateMask;

    if sum(cleanedAnalysisMask) >= 2
        cleanTimeMedian = median(chanSTD(cleanedAnalysisMask), 'omitnan');
        cleanTimeSigma = robustSigmaMAD(chanSTD(cleanedAnalysisMask));
        cleanTimeThresh = cleanTimeMedian + k * cleanTimeSigma;
    else
        cleanTimeMedian = globalTimeMedian;
        cleanTimeSigma = globalTimeSigma;
        cleanTimeThresh = globalTimeThresh;
    end

    %% -------------------- Lead-wise local thresholds --------------------
    leadThreshes = nan(1,numel(leads));

    for l = 1:numel(leads)
        leadCh = leads{l};
        leadCh = leadCh(leadCh >= 1 & leadCh <= nCh);

        leadGood = leadCh(analysisMask(leadCh));
        if numel(leadGood) < 2
            continue
        end

        leadSTD = chanSTD(leadGood);
        leadMedian = median(leadSTD, 'omitnan');
        leadSigma = robustSigmaMAD(leadSTD);

        if leadSigma > 3 * cleanTimeSigma
            leadThresh = cleanTimeMedian + 0.25 * leadSigma;
        else
            leadThresh = leadMedian + k * leadSigma;
        end

        leadThreshes(l) = leadThresh;

        leadCandidateMask(leadGood) = chanSTD(leadGood) > leadThresh;
    end

    candidateMask = (timeCandidateMask | spectralCandidateMask | leadCandidateMask) & analysisMask;

    %% -------------------- Spectral confirmation --------------------
    candidateMask = (timeCandidateMask | spectralCandidateMask | leadCandidateMask) & analysisMask;
    candidateIdx = find(candidateMask);
    goodMask = analysisMask & ~candidateMask;

    % Initialize outputs
    confirmedSpectralMask = false(1,nCh);

    if sum(goodMask) < 3 || isempty(candidateIdx) || isempty(freqVec)
        % Not enough reference → fallback
        confirmedSpectralMask = candidateMask;
    else
        goodEEG   = eegPower(goodMask);
        goodNoise = noiseFloor(goodMask);

        % --- Robust z-scores ---
        medEEG = median(goodEEG, 'omitnan');
        madEEG = 1.4826 * mad(goodEEG,1);
        medNoise = median(goodNoise, 'omitnan');
        madNoise = 1.4826 * mad(goodNoise,1);

        % Prevent division by zero
        madEEG(madEEG == 0) = eps;
        madNoise(madNoise == 0) = eps;
        zEEG   = abs((eegPower - medEEG) ./ madEEG);
        zNoise = abs((noiseFloor - medNoise) ./ madNoise);

        % --- MATLAB outlier detection (robust) ---
        outEEG   = isoutlier(eegPower,   'median', 'ThresholdFactor', k);
        outNoise = isoutlier(noiseFloor, 'median', 'ThresholdFactor', k);

        % Combine criteria
        confirmedSpectralMask = ((zEEG > k & zNoise > k)| outEEG | outNoise) & candidateMask;

    end
    %% -------------------- Final bad mask --------------------
    badMask = knownBadMask | impedanceBadMask | confirmedSpectralMask;

    badChannels = find(badMask);

    %% -------------------- Pack output info --------------------
    badChannelInfo = packInfo();

    %% -------------------- Plotting --------------------
    if graphicOut
        figure('Name', 'Bad Channel Diagnostics', 'Position', [300 200 1300 700]);

        tiledlayout(2,2, 'Padding','compact', 'TileSpacing','compact');

        nexttile
        plot(1:nCh, chanSTD, 'bo', 'DisplayName','Channel STD'); hold on
        yline(globalTimeMedian, '--k', 'DisplayName','Global median');
        yline(globalTimeThresh, '--r', 'DisplayName','Global threshold');
        for l = 1:numel(leads)
            if ~isnan(leadThreshes(l))
                plot(leads{l}, repmat(leadThreshes(l), size(leads{l})), '-.', ...
                    'DisplayName', sprintf('Lead %d threshold', l));
            end
        end
        scatter(find(badMask), chanSTD(badMask), 60, 'r', 'filled', 'DisplayName','Bad');
        xlabel('Channel');
        ylabel('STD');
        title('Time-domain channel STD');
        legend('Location','best');

        nexttile
        if ~isempty(freqVec)
            plot(1:nCh, spectralSTD, 'ko', 'DisplayName','Spectral STD'); hold on
            scatter(find(spectralCandidateMask), spectralSTD(spectralCandidateMask), 60, 'm', 'filled', ...
                'DisplayName','Spectral candidate');
            xlabel('Channel');
            ylabel('STD of log10 PSD');
            title('Spectral variability');
            legend('Location','best');
        else
            title('Spectral variability unavailable');
        end

        nexttile
        if ~isempty(freqVec)
            plot(1:nCh, eegPower, 'bo', 'DisplayName','EEG power'); hold on
            plot(1:nCh, noiseFloor, 'ro', 'DisplayName','Noise floor');
            scatter(find(badMask), eegPower(badMask), 60, 'k', 'filled', 'DisplayName','Bad');
            xlabel('Channel');
            ylabel('Mean log10 power');
            title('EEG power and noise floor');
            legend('Location','best');
        end

        nexttile
        if ~isempty(freqVec)
            hold on
            % Extract PSDs
            goodPSD = pxx(:, goodMask);
            badPSD  = pxx(:, badChannels);

            % Plot individual traces (faint)
            if ~isempty(goodPSD)
                plot(freqVec, goodPSD, 'Color', [0 0 1 0.15], 'HandleVisibility','off'); % light blue
            end

            if ~isempty(badPSD)
                plot(freqVec, badPSD, 'Color', [1 0 0 0.15],'HandleVisibility','off'); % light red
            end

            % Plot mean curves (bold)
            if ~isempty(goodPSD)
                plot(freqVec, mean(goodPSD,2,'omitnan'), ...
                    'b', 'LineWidth', 2, 'DisplayName','Good channels');
            end

            if ~isempty(badPSD)
                plot(freqVec, mean(badPSD,2,'omitnan'), ...
                    'r', 'LineWidth', 2, 'DisplayName','Bad channels');
            end
             % Log scale
            set(gca, 'YScale', 'log')
            xlabel('Frequency [Hz]')
            ylabel('Power')
            title('PSD comparison')
            legend('Location','best') 
        end

        fontsize(12, 'points');
    end

    %% -------------------- Nested packer --------------------
    function info = packInfo()
        info = struct();

        info.badMask = badMask;

        info.metrics = struct();
        info.metrics.channelSTD = chanSTD;
        if exist('spectralSTD','var'), info.metrics.spectralSTD = spectralSTD; else, info.metrics.spectralSTD = []; end
        if exist('eegPower','var'), info.metrics.eegPower = eegPower; else, info.metrics.eegPower = []; end
        if exist('noiseFloor','var'), info.metrics.noiseFloor = noiseFloor; else, info.metrics.noiseFloor = []; end

        info.thresholds = struct();
        info.thresholds.globalTimeMedian = globalTimeMedian;
        info.thresholds.globalTimeSigma = globalTimeSigma;
        info.thresholds.globalTimeThreshold = globalTimeThresh;
        if exist('cleanTimeMedian','var'), info.thresholds.cleanTimeMedian = cleanTimeMedian; else, info.thresholds.cleanTimeMedian = []; end
        if exist('cleanTimeSigma','var'), info.thresholds.cleanTimeSigma = cleanTimeSigma; else, info.thresholds.cleanTimeSigma = []; end
        if exist('cleanTimeThresh','var'), info.thresholds.cleanTimeThreshold = cleanTimeThresh; else, info.thresholds.cleanTimeThreshold = []; end
        info.thresholds.leadThresholds = leadThreshes;

        thresholds = struct();
        thresholds.globalTimeThreshold = globalTimeThresh;
        thresholds.cleanTimeThreshold = cleanTimeThresh;
        thresholds.leadThresholds = leadThreshes;
        thresholds.k = k;

        info.parameters = struct();
        info.parameters.k = k;
        info.parameters.FS = FS;
        info.parameters.construct = construct;
        info.parameters.leads = leads;
        info.parameters.impedanceThreshold = 10000;
        info.parameters.eegBandHz = [1 min(150, nyq)];
        info.parameters.noiseBandHz = [noiseLow nyq];
    end
end


%% Helper functions

function sig = robustSigmaMAD(x)
    x = x(:);
    x = x(~isnan(x) & isfinite(x));

    if isempty(x)
        sig = NaN;
        return
    end

    medx = median(x);
    madx = median(abs(x - medx));

    sig = 1.4826 * madx;

    if sig == 0 || isnan(sig)
        sig = std(x, 0, 'omitnan');
    end

    if sig == 0 || isnan(sig)
        sig = eps;
    end
end


function chIdx = parseChannelIndex(x)
    chIdx = [];

    if isnumeric(x) && isscalar(x)
        chIdx = double(x);
        return
    end

    if iscell(x) && ~isempty(x)
        x = x{1};
    end

    if ischar(x) || isstring(x)
        tok = regexp(char(string(x)), '\d+', 'match', 'once');
        if ~isempty(tok)
            chIdx = str2double(tok);
        end
    end
end

