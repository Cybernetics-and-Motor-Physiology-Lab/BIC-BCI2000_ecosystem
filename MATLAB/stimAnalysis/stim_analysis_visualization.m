clear all; close all; clc
%% Add Paths
addpath(genpath('./functions'));
addpath(genpath('./colormaps'));
addpath('.../violinplot'); % https://github.com/bastibe/Violinplot-Matlab.git

%% Load KJM colormaps
rwb = load('rwb_colormap.mat');
rgb = load('RYGyCB.mat');

%% Channel map ( will be saved in next iteration in data structure)
channelLeadMap = dictionary( ...
    'displaced', struct('indices', 1:8, 'labels',  "Ch" + (1:8)), ...
    'occipital',struct('indices', 9:16, 'labels',  "Ch" + (9:16)), ...
    'thalamus',struct('indices', 17:24, 'labels',  "Ch" + (17:24)),...
    'lobus-piriformis', struct('indices', 25:32, 'labels',  "Ch" + (25:32)) );
%% Load the processed data
[fname, fpth] = uigetfile(['/path/to/data' filesep '*.mat'], ...
    'Select the INPUT .mat FILE','MultiSelect','off');      % UI selectfile to open
load([fpth filesep fname]);

%% Extratc main parameters
fld_names = fieldnames(data);

%% Iterate over the files and plot average stim spectogram and PSD
% speciy channel to plot
ch = 21;

% Deine spectral parameters
freqSpec = 1:500;
win = 5;
noverlap = 0.5;

for f = 1:numel(fld_names)
    if ~isfield(data.(fld_names{f}),'StimInfo') % Skip baseline
        continue
    end
    stimInfo = data.(fld_names{f}).StimInfo;
    stimPulse = stimInfo.StimPulses;
    stimT = (5*stimPulse{3,1} + 2*stimPulse{4,1} + stimPulse{5,1})/1e6; % Convert from µs to s
    stimF = 1/stimT;    % convert to f

    stimSeqMask= stimInfo.StimSequence.Sequence~=0; % Find non-zero sequences (stim)
    stimOnsets = stimInfo.StimSequence.SequenceOnset(stimSeqMask);
    stimDurs = stimInfo.StimSequence.SequenceDuration(stimSeqMask);

    FS = data.(fld_names{f}).SamplingRate;
    sig = data.(fld_names{f}).Signal;

    figure('Position', [200, 100, 1000, 800]);
    hold on
    stimPSDs = zeros(numel(freqSpec), sum(stimSeqMask));
    for s = 1:sum(stimSeqMask)
        stimIdx = stimOnsets(s):stimOnsets(s)+stimDurs-1;
        [stimPSDs(:,s),~] = pwelch(sig(stimIdx,ch),win*FS,noverlap*(win*FS),freqSpec,FS);
        plot(freqSpec,stimPSDs(:,s), 'Color', [0 0 0 0.8], 'LineWidth', 1.5, ...
            'DisplayName', sprintf('StimTrial%d', s));
    end
    plot(freqSpec,mean(stimPSDs,2), 'Color', [1 0 0], 'DisplayName', 'Average', 'LineWidth', 2.5);

    title(sprintf('Stim Frequnecy: %d Hz', stimF));
    sgtitle(sprintf('File: %s',fld_names{f}));
    xlim tight;
    xlabel('Frequency (Hz)');
    legend(legendUnq());
    grid on;
    fontsize(16,'points');
    yscale log
end

clearvars -except fld_names data rgb rwb fname fpth channelLeadMap

%% Iterate over the files and plot number of spikes
 % Spike counts and normalized counts
spikeCounts = nan(1,numel(fld_names));
norm_nSpike = nan(1,numel(fld_names));
totalSpikePerChannel = struct();

for f = 1:numel(fld_names)
    spikeMap = containers.Map('KeyType', 'char', 'ValueType', 'double');  % aggregate spike counts
    FS = data.(fld_names{f}).SamplingRate;
    spikeDetect = data.(fld_names{f}).SpikeDetection;
    badChnls = data.(fld_names{f}).BadChannels.BadChannels;
    nCh = 32 - numel(badChnls);
    chanNames = spikeDetect.Montage.ChannelNames;      % Retrieve the channel names
    spikeChannels = spikeDetect.Channelinformation(spikeDetect.RFOMP.feature.Pred==1); % Retrieve channel with spikes
    % Remove spikes detected on bad channels
    spikeChannels(ismember(spikeChannels,badChnls)) = [];
    nSpike = numel(spikeChannels);

    % Accumulate into map
    for c = 1:numel(spikeChannels)
        ch = chanNames{spikeChannels(c)};
        if spikeMap.isKey(ch)
            spikeMap(ch) = spikeMap(ch) + 1;
        else
            spikeMap(ch) = 1;
        end
    end

    if isfield(data.(fld_names{f}), 'StimInfo') % Stim Recordings
        seqC = data.(fld_names{f}).StimInfo.StimSequence.Sequence;
        seqDur = data.(fld_names{f}).StimInfo.StimSequence.SequenceDuration;
        dur = minutes(seconds(sum(seqDur(seqC == 0))/FS)); % Duration of the non-stim periods
    else % Baseline
        dur = minutes(duration(data.(fld_names{f}).Duration,'Format','hh:mm:ss'));
    end

    spikeCounts(f) = nSpike;
    norm_nSpike(f) = nSpike / (dur*nCh); % Normalized spike count by dividing through number of channels and duration of recording

    % Convert map to table
    allKeys = keys(spikeMap);
    allValues = values(spikeMap);
    spikeTbl = table(allKeys', cell2mat(allValues'), ...
        'VariableNames', {'Channel', 'SpikeCount'});
    totalSpikePerChannel.(fld_names{f}).SpikeTable = spikeTbl;

end

figure('Position', [200, 100, 1200, 800]);
bar(norm_nSpike)
title('Spike rate across the recordings')
ylabel('Normalized spike rate [nSpikes / (nCh*minutes)]')
xticklabels(fld_names)
ax = gca;
ax.TickLabelInterpreter = 'none';

clearvars -except fld_names data rgb rwb fname fpth norm_nSpike totalSpikePerChannel channelLeadMap

%% Plot the heatmaps of temporal features
temp_features = {'CrossCorrelation'};
cmap = 'hot';

for f = 1:numel(fld_names)
    ch_names = data.(fld_names{f}).ChannelNames;
    allTrials = fieldnames(data.(fld_names{f}).TemporalFeatures);
    isTrialField = ~cellfun(@isempty, regexp(allTrials, '^Trial\d+$', 'once')); % Keep only fields named Trial%d
    trials = allTrials(isTrialField);
    trialNums = cellfun(@(x) sscanf(x, 'Trial%d'), trials); % Sort trials numerically
    [trialNums, sortIdx] = sort(trialNums);
    trials = trials(sortIdx);
    nTrials = numel(trials);

    for ft = 1:numel(temp_features)
        current_feature = temp_features{ft};
        figure('Position', [200, 100, 1400, 1000]);
        % Layout definition
        if nTrials > 2
            nLeftRows = ceil(sqrt(nTrials));
            nLeftCols = ceil(nTrials / nLeftRows);
            % Total layout: left side trial plots, right side average plot
            tl = tiledlayout(nLeftRows, nLeftCols + 1, 'TileSpacing', 'compact', ...
                'Padding', 'compact');
            avgTileIdx = nLeftCols + 1;
        elseif nTrials == 2
            nLeftRows = 2;
            nLeftCols = 1;
            tl = tiledlayout(2, 2, 'TileSpacing', 'compact','Padding', 'compact');
            avgTileIdx = 2;
        else
            tl = tiledlayout(1, 1, 'TileSpacing', 'compact', 'Padding', 'compact');
        end

        sgtitle(sprintf('%s - %s', fld_names{f}, current_feature), ...
            'Interpreter', 'none');

        % Collect matrices for averaging
        featureMats = [];

        for tr = 1:nTrials
            currTrial = trials{tr};
            if ~isfield(data.(fld_names{f}).TemporalFeatures.(currTrial), current_feature)
                warning('Missing %s in %s.', current_feature, currTrial);
                continue
            end
            M = data.(fld_names{f}).TemporalFeatures.(currTrial).(current_feature);
            if isempty(M)
                continue
            end
            featureMats(:,:,end+1) = M; %#ok<SAGROW>
            % Select tile
            if nTrials > 2
                rowIdx = ceil(tr / nLeftCols);
                colIdx = mod(tr-1, nLeftCols) + 1;
                tileIdx = (rowIdx-1) * (nLeftCols + 1) + colIdx;
                nexttile(tl, tileIdx);
            elseif nTrials == 2
                nexttile(tl, (tr-1)*2 + 1);
            else
                nexttile(tl, 1);
            end

            imagesc(M)
            axis square
            colormap(cmap)
            colorbar
            title(sprintf('Trial %d', trialNums(tr)))
            % xticks(1:numel(ch_names))
            % yticks(1:numel(ch_names))
            % xticklabels(ch_names)
            % yticklabels(ch_names)
        end

        % Average plot
        if nTrials > 1 && ~isempty(featureMats)
            avgM = mean(featureMats, 3, 'omitnan');
            if nTrials > 2
                % Right side spans all rows
                nexttile(tl, avgTileIdx, [nLeftRows 1]);
            elseif nTrials == 2
                nexttile(tl, 2, [2 1]);
            end
            imagesc(avgM)
            axis square
            colormap(cmap)
            colorbar
            title('Average')
            % xticks(1:numel(ch_names))
            % yticks(1:numel(ch_names))
            % xticklabels(ch_names)
            % yticklabels(ch_names)
        end
        fontsize(16,'points');
    end
end

clearvars -except fld_names data rgb rwb fname fpth norm_nSpike totalSpikePerChannel channelLeadMap

%% Plot Spectral Features
spect_features = {'BandPower', 'PSD'}; % 'BandCoherence'
hm_cmap = hot;

runSubfieldAvg = struct();
leadMap = keys(channelLeadMap);
leadColors = lines(numel(leadMap));

for f = 1:numel(fld_names)

    fname = fld_names{f};
    ch_names = data.(fname).ChannelNames;

    allTrials = fieldnames(data.(fname).SpectralFeatures);
    isTrialField = ~cellfun(@isempty, regexp(allTrials, '^Trial\d+$', 'once'));
    trials = allTrials(isTrialField);

    trialNums = cellfun(@(x) sscanf(x, 'Trial%d'), trials);
    [trialNums, sortIdx] = sort(trialNums);
    trials = trials(sortIdx);
    nTrials = numel(trials);

    % Good channels
    badChnls = data.(fname).BadChannels.BadChannels;
    goodChnls = 1:numel(ch_names);
    goodChnls(badChnls) = [];

    for ft = 1:numel(spect_features)
        current_feature = spect_features{ft};

        % Get first available trial to inspect feature structure
        firstTrial = trials{1};

        if ~isfield(data.(fname).SpectralFeatures.(firstTrial), current_feature)
            warning('Missing feature %s in %s.', current_feature, fname);
            continue
        end

        featureObj = data.(fname).SpectralFeatures.(firstTrial).(current_feature);

        hasSubfields = isstruct(featureObj);

        featureMats = [];
        trialSubfieldAvg = struct();

        % Create figure and layout
        figure('Position', [200, 100, 1500, 950]);

        if hasSubfields
            subfields = fieldnames(featureObj);
            nSubfields = numel(subfields);
            tl = createSubfieldLayout(nSubfields, nTrials);
        else
            tl = createTrialLayout(nTrials);
        end

        sgtitle(sprintf('%s | %s', fname, current_feature), ...
            'Interpreter', 'none', ...
            'FontSize', 18, ...
            'FontWeight', 'bold');

        % Plot individual trials
        for tr = 1:nTrials
            currTrial = trials{tr};
            if ~isfield(data.(fname).SpectralFeatures.(currTrial), current_feature)
                warning('Missing %s in %s.', current_feature, currTrial);
                continue
            end

            switch current_feature
                case 'BandPower'
                    for b = 1:nSubfields
                        currBand = subfields{b};
                        bandData = data.(fname).SpectralFeatures.(currTrial). ...
                            (current_feature).(currBand);   % expected [time x channels]

                        ax = nexttile(tl, tileIndexSubfield(b, tr, nTrials));
                        hold(ax, 'on');

                        for l = 1:numel(leadMap)
                            leadIdxs = channelLeadMap(leadMap(l)).indices;
                            leadIdxs = intersect(leadIdxs, goodChnls);

                            if isempty(leadIdxs)
                                continue
                            end

                            thisColor = leadColors(l,:);

                            % Individual channels
                            plot(ax, bandData(:,leadIdxs), ...
                                'Color', [thisColor 0.35], ...
                                'LineWidth', 0.75, ...
                                'HandleVisibility', 'off');

                            % Lead average
                            plot(ax, mean(bandData(:,leadIdxs), 2, 'omitnan'), ...
                                'Color', thisColor, ...
                                'LineWidth', 2, ...
                                'DisplayName', string(leadMap(l)));
                        end

                        title(ax, sprintf('Trial %d | %s', trialNums(tr), currBand), ...
                            'Interpreter', 'none');
                        ylabel(ax, 'Power [\muV^2]');
                        xlabel(ax, 'Time bin');
                        grid(ax, 'on');

                        if tr == 1 && b == 1
                            legend(ax, 'Location', 'best');
                        end

                        trialSubfieldAvg.(currBand)(tr,1) = ...
                            mean(bandData(:,goodChnls), 'all', 'omitnan');
                    end

                case 'BandCoherence'
                    for b = 1:nSubfields

                        currBand = subfields{b};
                        bandData = data.(fname).SpectralFeatures.(currTrial). ...
                            (current_feature).(currBand);   % expected [channels x channels]

                        ax = nexttile(tl, tileIndexSubfield(b, tr, nTrials));

                        imagesc(ax, bandData);
                        axis(ax, 'square');
                        colormap(ax, hm_cmap);
                        colorbar(ax);

                        title(ax, sprintf('Trial %d | %s', trialNums(tr), currBand), ...
                            'Interpreter', 'none');

                        xlabel(ax, 'Channel');
                        ylabel(ax, 'Channel');

                        xticks(ax, 1:numel(ch_names));
                        yticks(ax, 1:numel(ch_names));
                        xticklabels(ax, ch_names);
                        yticklabels(ax, ch_names);
                        ax.TickLabelInterpreter = 'none';
                        ax.FontSize = 8;

                        trialSubfieldAvg.(currBand)(:,:,tr) = bandData;
                    end

                case 'PSD'
                    psd = data.(fname).SpectralFeatures.(currTrial).(current_feature);

                    freqBins = getFreqBins(data.(fname).SpectralFeatures.(currTrial));

                    ax = nexttile(tl, tr);
                    hold(ax, 'on');

                    for l = 1:numel(leadMap)
                        leadIdxs = channelLeadMap(leadMap(l)).indices;
                        leadIdxs = intersect(leadIdxs, goodChnls);

                        if isempty(leadIdxs)
                            continue
                        end

                        thisColor = leadColors(l,:);

                        % Individual channels
                        plot(ax, freqBins, psd(:,leadIdxs), ...
                            'Color', [thisColor 0.25], ...
                            'LineWidth', 0.75, ...
                            'HandleVisibility', 'off');

                        % Lead average
                        plot(ax, freqBins, mean(psd(:,leadIdxs), 2, 'omitnan'), ...
                            'Color', thisColor, ...
                            'LineWidth', 2, ...
                            'DisplayName', string(leadMap(l)));
                    end

                    set(ax, 'YScale', 'log');
                    ylabel(ax, 'PSD [\muV^2/Hz]');
                    xlabel(ax, 'Frequency [Hz]');
                    title(ax, sprintf('Trial %d PSD', trialNums(tr)));
                    xlim(ax, [min(freqBins) max(freqBins)]);
                    grid(ax, 'on');

                    if tr == 1
                        legend(ax, 'Location', 'best');
                    end

                    featureMats(:,:,tr) = psd; %#ok<SAGROW>
            end
        end

        %% Add average tile for PSD if available
        if strcmp(current_feature, 'PSD') && ~isempty(featureMats) && nTrials > 1

            avgTileIdx = getAverageTileIndex(nTrials);
            ax = nexttile(tl, avgTileIdx);
            hold(ax, 'on');

            avgPSD = mean(featureMats, 3, 'omitnan');
            freqBins = getFreqBins(data.(fname).SpectralFeatures.(trials{1}));

            for l = 1:numel(leadMap)
                leadIdxs = channelLeadMap(leadMap(l)).indices;
                leadIdxs = intersect(leadIdxs, goodChnls);

                if isempty(leadIdxs)
                    continue
                end

                thisColor = leadColors(l,:);

                plot(ax, freqBins, mean(avgPSD(:,leadIdxs), 2, 'omitnan'), ...
                    'Color', thisColor, ...
                    'LineWidth', 2.5, ...
                    'DisplayName', string(leadMap(l)));
            end

            set(ax, 'YScale', 'log');
            title(ax, 'Average PSD across trials');
            xlabel(ax, 'Frequency [Hz]');
            ylabel(ax, 'PSD [\muV^2/Hz]');
            grid(ax, 'on');
            legend(ax, 'Location', 'best');
        end

        %% Save values for group-level plots
        if hasSubfields
            runSubfieldAvg.(fname).(current_feature) = trialSubfieldAvg;
        else
            runSubfieldAvg.(fname).(current_feature) = featureMats;
        end
    end
end

% Plot averaged figures across runs/files
for ft = 1:numel(spect_features)

    current_feature = spect_features{ft};

    switch current_feature

        case 'BandPower'
            plotBandPowerRunSummary(runSubfieldAvg, fld_names, current_feature);

        case 'BandCoherence'
            plotBandCoherenceRunSummary(runSubfieldAvg, fld_names, current_feature, hm_cmap, ch_names);

        case 'PSD'
            plotPSDRunSummary(runSubfieldAvg, fld_names, current_feature, channelLeadMap);
    end
end
function tl = createTrialLayout(nTrials)
    if nTrials > 2
        nLeftRows = ceil(sqrt(nTrials));
        nLeftCols = ceil(nTrials / nLeftRows);

        tl = tiledlayout(nLeftRows, nLeftCols + 1, ...
            'TileSpacing', 'compact', ...
            'Padding', 'compact');

    elseif nTrials == 2
        tl = tiledlayout(2, 2, ...
            'TileSpacing', 'compact', ...
            'Padding', 'compact');

    else
        tl = tiledlayout(1, 1, ...
            'TileSpacing', 'compact', ...
            'Padding', 'compact');
    end
end


function tl = createSubfieldLayout(nSubfields, nTrials)
    tl = tiledlayout(nSubfields, nTrials, ...
        'TileSpacing', 'compact', ...
        'Padding', 'compact');
end


function idx = tileIndexSubfield(subfieldIdx, trialIdx, nTrials)
    idx = (subfieldIdx - 1) * nTrials + trialIdx;
end


function avgTileIdx = getAverageTileIndex(nTrials)
    if nTrials > 2
        nLeftRows = ceil(sqrt(nTrials));
        nLeftCols = ceil(nTrials / nLeftRows);
        avgTileIdx = nLeftCols + 1;
    elseif nTrials == 2
        avgTileIdx = 2;
    else
        avgTileIdx = 1;
    end
end


function freqBins = getFreqBins(trialStruct)
    if isfield(trialStruct, 'info') && isfield(trialStruct.info, 'Frequencies')
        freqBins = trialStruct.info.Frequencies;
    elseif isfield(trialStruct, 'Info') && isfield(trialStruct.Info, 'Frequencies')
        freqBins = trialStruct.Info.Frequencies;
    else
        error('Could not find frequency vector in trialStruct.info.Frequencies.');
    end

    freqBins = freqBins(:);
end


function plotBandPowerRunSummary(runSubfieldAvg, fld_names, featureName)

    firstRun = fld_names{1};
    bandNames = fieldnames(runSubfieldAvg.(firstRun).(featureName));

    figure('Position', [200, 100, 1300, 700]);
    tl = tiledlayout(1, numel(bandNames), ...
        'TileSpacing', 'compact', ...
        'Padding', 'compact');

    sgtitle('BandPower summary across runs', ...
        'FontSize', 18, ...
        'FontWeight', 'bold');

    for b = 1:numel(bandNames)
        bandName = bandNames{b};
        ax = nexttile(tl);
        hold(ax, 'on');

        allVals = [];
        cats = {};

        for f = 1:numel(fld_names)
            fname = fld_names{f};

            if ~isfield(runSubfieldAvg.(fname).(featureName), bandName)
                continue
            end

            vals = runSubfieldAvg.(fname).(featureName).(bandName);
            vals = vals(:);

            allVals = [allVals; vals]; %#ok<AGROW>
            cats = [cats; repmat({fname}, numel(vals), 1)]; %#ok<AGROW>
        end

        boxchart(ax, categorical(cats), allVals);
        title(ax, bandName, 'Interpreter', 'none');
        ylabel(ax, 'Mean band power [\muV^2]');
        ax.TickLabelInterpreter = 'none';
        grid(ax, 'on');
    end
end


function plotBandCoherenceRunSummary(runSubfieldAvg, fld_names, featureName, hm_cmap, ch_names)

    firstRun = fld_names{1};
    bandNames = fieldnames(runSubfieldAvg.(firstRun).(featureName));

    for b = 1:numel(bandNames)
        bandName = bandNames{b};

        figure('Position', [200, 100, 1400, 800]);
        tl = tiledlayout(1, numel(fld_names), ...
            'TileSpacing', 'compact', ...
            'Padding', 'compact');

        sgtitle(sprintf('Average BandCoherence: %s', bandName), ...
            'Interpreter', 'none', ...
            'FontSize', 18, ...
            'FontWeight', 'bold');

        for f = 1:numel(fld_names)
            fname = fld_names{f};

            ax = nexttile(tl);

            if ~isfield(runSubfieldAvg.(fname).(featureName), bandName)
                title(ax, sprintf('%s missing', fname), 'Interpreter', 'none');
                continue
            end

            M = mean(runSubfieldAvg.(fname).(featureName).(bandName), 3, 'omitnan');

            imagesc(ax, M);
            axis(ax, 'square');
            colormap(ax, hm_cmap);
            colorbar(ax);

            title(ax, fname, 'Interpreter', 'none');
            xlabel(ax, 'Channel');
            ylabel(ax, 'Channel');

            xticks(ax, 1:numel(ch_names));
            yticks(ax, 1:numel(ch_names));
            xticklabels(ax, ch_names);
            yticklabels(ax, ch_names);
            ax.TickLabelInterpreter = 'none';
            ax.FontSize = 8;
        end
    end
end


function plotPSDRunSummary(runSubfieldAvg, fld_names, featureName, channelLeadMap)

    leadMap = keys(channelLeadMap);
    leadColors = lines(numel(leadMap));

    figure('Position', [200, 100, 1200, 800]);
    ax = axes;
    hold(ax, 'on');

    for f = 1:numel(fld_names)
        fname = fld_names{f};

        if ~isfield(runSubfieldAvg.(fname), featureName)
            continue
        end

        psdStack = runSubfieldAvg.(fname).(featureName);

        if isempty(psdStack)
            continue
        end

        avgPSD = mean(psdStack, 3, 'omitnan');

        % Frequency axis is unavailable from runSubfieldAvg alone,
        % so use bin number here unless you pass freqBins into this helper.
        freqBins = 1:size(avgPSD,1);

        for l = 1:numel(leadMap)
            leadIdxs = channelLeadMap(leadMap(l)).indices;
            leadIdxs = leadIdxs(leadIdxs <= size(avgPSD,2));

            thisColor = leadColors(l,:);

            plot(ax, freqBins, mean(avgPSD(:,leadIdxs), 2, 'omitnan'), ...
                'Color', thisColor, ...
                'LineWidth', 1.5, ...
                'DisplayName', sprintf('%s | %s', fname, string(leadMap(l))));
        end
    end

    set(ax, 'YScale', 'log');
    xlabel(ax, 'Frequency bin');
    ylabel(ax, 'PSD [\muV^2/Hz]');
    title(ax, 'Average PSD across trials and runs');
    grid(ax, 'on');
    legend(ax, 'Location', 'bestoutside', 'Interpreter', 'none');
end