function h = plotTrialSpectralBrowser(signal, trialSpectStruct, fs, channelNames)
% plotTrialSpectralBrowser Interactive browser for signal, TF map, and PIB.
% Inputs
%   signal           - [samples x channels] signal matrix
%   trialSpectStruct - struct with fields:
%                        TFmap [freq x time x channel]
%                        TFtime
%                        freqBins
%                        trialVec
%                        trialStimCode
%                        validTrial
%                        trialStartIdx
%                        trialStopIdx
%   fs               - sampling rate in Hz
%   channelNames     - optional cell/string array of channel names
%
% Controls
%   Left / Right arrows     - change channel

    %% Validate inputs
    if nargin < 4 || isempty(channelNames)
        channelNames = arrayfun(@(ch) sprintf('Ch %d', ch), ...
            1:size(signal,2), 'UniformOutput', false);
    end

    requiredFields = {'TFmap','TFtime','freqBins','trialVec', ...
        'trialStimCode','validTrial','trialStartIdx','trialStopIdx'};

    if ~isstruct(trialSpectStruct)
        error('trialSpectStruct must be a nonempty struct.');
    end

    for i = 1:numel(requiredFields)
        if ~isfield(trialSpectStruct, requiredFields{i})
            error('trialSpectStruct is missing required field "%s".', requiredFields{i});
        end
    end

    signal = double(signal);
    channelNames = cellstr(string(channelNames));

    nSamples = size(signal,1);
    nCh = size(signal,2);

    if numel(channelNames) ~= nCh
        error('Number of channelNames must match number of signal channels.');
    end

    TFmap = trialSpectStruct.TFmap;
    TFtime = trialSpectStruct.TFtime(:);
    freqBins = trialSpectStruct.freqBins(:);
    meanPSD = trialSpectStruct.meanPSD;

    if size(TFmap,3) ~= nCh
        error('TFmap third dimension must match number of signal channels.');
    end

    %% State
    tSignal = (0:nSamples-1) ./ fs;
    maxTime = tSignal(end);

    state.ch = 1;
    state.fMin = 15;
    state.fMax = 30;

    %% GUI
    fig = uifigure('Name', 'Trial Spectral Browser', ...
        'Position', [100 100 1300 850], ...
        'KeyPressFcn', @keyPressHandler);

    mainGL = uigridlayout(fig, [2 1]);
    mainGL.RowHeight = {48, '1x'};
    mainGL.Padding = [10 10 10 10];
    mainGL.RowSpacing = 8;

    %% Top controls
    topControls = uigridlayout(mainGL, [1 9]);
    topControls.ColumnWidth = {70, 70, 95, 80, 95, 80, 120, 90, '1x'};
    topControls.Padding = [0 0 0 0];
    topControls.ColumnSpacing = 8;

    uilabel(topControls, 'Text', 'Channel');

    chEdit = uieditfield(topControls, 'numeric', ...
        'Value', state.ch, ...
        'Limits', [1 nCh], ...
        'RoundFractionalValues', true, ...
        'ValueChangedFcn', @channelEditCallback);

    uilabel(topControls, 'Text', 'PIB min [Hz]');

    fMinField = uieditfield(topControls, 'numeric', ...
        'Value', state.fMin, ...
        'Limits', [min(freqBins) max(freqBins)], ...
        'ValueChangedFcn', @freqFieldCallback);

    uilabel(topControls, 'Text', 'PIB max [Hz]');

    fMaxField = uieditfield(topControls, 'numeric', ...
        'Value', state.fMax, ...
        'Limits', [min(freqBins) max(freqBins)], ...
        'ValueChangedFcn', @freqFieldCallback);

    normalizeTFBox = uicheckbox(topControls, ...
        'Text', 'Normalize TF', ...
        'Value', false, ...
        'ValueChangedFcn', @(~,~) updatePlot());

    exportFigBtn = uibutton(topControls, ...
        'Text', 'Open Figure', ...
        'ButtonPushedFcn', @(~,~) openCurrentAxesInFigure());

    %% Plot panel containing only axes
    plotPanel = uipanel(mainGL, ...
        'BorderType', 'none', ...
        'BackgroundColor', [1 1 1]);

    plotGL = uigridlayout(plotPanel, [3 1]);
    plotGL.RowHeight = {'1x', '1x', '1x'};
    plotGL.Padding = [10 10 10 10];
    plotGL.RowSpacing = 10;

    axSig = uiaxes(plotGL);
    axTF  = uiaxes(plotGL);
    axPIB = uiaxes(plotGL);

    %% Initial plot
    updatePlot();

    %% Output handles
    h = struct();
    h.fig = fig;
    h.plotPanel = plotPanel;
    h.axSignal = axSig;
    h.axTF = axTF;
    h.axPIB = axPIB;
    h.channelEdit = chEdit;
    h.fMinField = fMinField;
    h.fMaxField = fMaxField;
    h.normalizeTFBox = normalizeTFBox;

    %% Callback functions
    function keyPressHandler(~, event)
        switch event.Key
            case 'rightarrow'
                state.ch = max(1, state.ch + 1);
            case 'leftarrow'
                state.ch = min(nCh, state.ch - 1);
        end

        chEdit.Value = state.ch;
        updatePlot();
    end

    function channelEditCallback(src, ~)
        state.ch = round(src.Value);
        updatePlot();
    end

    function freqFieldCallback(~, ~)
        state.fMin = fMinField.Value;
        state.fMax = fMaxField.Value;

        if state.fMin >= state.fMax
            uialert(fig, 'PIB min must be smaller than PIB max.', ...
                'Invalid frequency band');

            state.fMin = 15;
            state.fMax = 30;

            fMinField.Value = state.fMin;
            fMaxField.Value = state.fMax;
        end

        updatePlot();
    end

    function openCurrentAxesInFigure()
        ch = state.ch;
        exportFig = figure( ...
            'Name', sprintf('Trial Spectral Browser - %s', channelNames{ch}), ...
            'Color', 'w', ...
            'Position', [100 100 1300 850]);
        tl = tiledlayout(exportFig, 3, 1, ...
            'TileSpacing', 'compact', ...
            'Padding', 'compact');

        % Copy signal axis
        axNew1 = nexttile(tl);
        copyAxesContent(axSig, axNew1);
        title(axNew1, axSig.Title.String, 'Interpreter', 'none');
        ylabel(axNew1, axSig.YLabel.String);
        grid(axNew1, 'on');

        % Copy TF axis
        axNew2 = nexttile(tl);
        copyAxesContent(axTF, axNew2);
        title(axNew2, axTF.Title.String, 'Interpreter', 'none');
        ylabel(axNew2, axTF.YLabel.String);
        colormap(axNew2, turbo);
        colorbar(axNew2);
        set(axNew2, 'YDir', 'normal');

        % Copy PIB axis
        axNew3 = nexttile(tl);
        copyAxesContent(axPIB, axNew3);
        title(axNew3, axPIB.Title.String, 'Interpreter', 'none');
        ylabel(axNew3, axPIB.YLabel.String);
        xlabel(axNew3, axPIB.XLabel.String);
        grid(axNew3, 'on');

        linkaxes([axNew1 axNew2 axNew3], 'x');
    end

    %% Main plotting function
    function updatePlot()
        ch = state.ch;
        t0 = 0;
        t1 = maxTime;

        %% Signal plot
        cla(axSig, 'reset')
        hold(axSig, 'on')

        sigFull = signal(:,ch);
        % sigYLims = prctile(sigFull(isfinite(sigFull)), [1 99]);
        [sigYLims(1) sigYLims(2)] = bounds(sigFull);
        if sigYLims(1) == sigYLims(2)
            sigYLims = sigYLims + [-1 1];
        end

        plot(axSig, tSignal, sigFull, ...
            'Color', [0.1 0.1 0.1], ...
            'LineWidth', 1);

        ylim(axSig, sigYLims)
        patchTrials(axSig, t0, t1, sigYLims)

        ylabel(axSig, 'Amplitude [µV]')
        title(axSig, sprintf('Signal: %s', channelNames{ch}), ...
            'Interpreter', 'none')
        grid(axSig, 'on')
        xlim(axSig, [t0 t1])

        %% TF map
        cla(axTF, 'reset')

        tfData = TFmap(:,:,ch);

        if normalizeTFBox.Value
            tfPlot = pow2db(tfData + eps) - pow2db(meanPSD(:,ch) + eps);
            tfTitle = 'Normalized time-frequency map';
            colorLimits = [-10 10];
        else
            tfPlot = pow2db(tfData + eps);
            tfTitle = 'Time-frequency map [dB PSD]';
            colorLimits = [-30 30];
        end

        imagesc(axTF, TFtime, freqBins, tfPlot)
        set(axTF, 'YDir', 'normal')
        clim(axTF, colorLimits)
        colormap(axTF, turbo)
        cBar = colorbar(axTF);
        cBar.Label.String = 'Power [dB]';
        ylabel(axTF, 'Frequency [Hz]')
        title(axTF, tfTitle)
        xlim(axTF, [t0 t1])
        ylim(axTF, [min(freqBins) max(freqBins)])

        hold(axTF, 'on')
        plotTrialLines(axTF, t0, t1)

        %% PIB plot
        cla(axPIB, 'reset')
        hold(axPIB, 'on')

        fIdx = freqBins >= state.fMin & freqBins < state.fMax;

        if any(fIdx)
            pib = squeeze(trapz(freqBins(fIdx), TFmap(fIdx,:,ch), 1));
            pib = pib(:);

            pibYLims = prctile(pib(isfinite(pib)), [1 99]);

            if pibYLims(1) == pibYLims(2)
                pibYLims = pibYLims + [-1 1];
            end

            plot(axPIB, TFtime, pib, ...
                'LineWidth', 1.5, ...
                'Color', [0 0.4470 0.7410]);

            ylim(axPIB, pibYLims)
            patchTrials(axPIB, t0, t1, pibYLims)

            ylabel(axPIB, 'Power [µV^(2)]')
            title(axPIB, sprintf('Power in band %d-%d Hz', ...
                state.fMin, state.fMax))
            grid(axPIB, 'on')
            xlim(axPIB, [t0 t1])
        else
            title(axPIB, 'No frequency bins in selected band')
        end

        xlabel(axPIB, 'Time [s]')

        linkaxes([axSig axTF axPIB], 'x')

        % Match signal and PIB axes width to TF axis, which has a colorbar.
        drawnow
        tfPos = axTF.Position;
        sigPos = axSig.Position;
        pibPos = axPIB.Position;
        axSig.Position = [tfPos(1), sigPos(2), tfPos(3), sigPos(4)];
        axPIB.Position = [tfPos(1), pibPos(2), tfPos(3), pibPos(4)];

        %% Create legend
        hold(axPIB, 'on')

        legendColors = {
            [0.5 0.5 0.5],          'StimCode 0';
            [0 0.4470 0.7410],      'StimCode 1';
            [0.9290 0.6940 0.1250], 'StimCode 2';
            [0.4660 0.6740 0.1880], 'StimCode 3';
            [1      0.5333 0],      'StimCode 4';
            [0.4940 0.1840 0.5560], 'StimCode 5';
            [1      0.5    1],      'Other StimCode'
            [1 0 0],                'Invalid Trial'
            };

        lgdHandles = gobjects(size(legendColors,1),1);

        for k = 1:size(legendColors,1)
            lgdHandles(k) = patch(axPIB, ...
                nan, nan, legendColors{k,1}, ...
                'FaceAlpha', 0.7, ...
                'EdgeColor', 'none', ...
                'Visible', 'off');
        end

        stimLegend = legend(axPIB, lgdHandles, legendColors(:,2), ...
            'Orientation', 'horizontal', ...
            'Location', 'southoutside');

    end


%% Trial patch helper
    function patchTrials(ax, t0, t1, yl)
        delete(findall(ax, 'Type', 'patch'))

        for tr = 1:numel(trialSpectStruct.trialStimCode)
            trStart = trialSpectStruct.trialStartIdx(tr) / fs;
            trStop  = trialSpectStruct.trialStopIdx(tr) / fs;

            if trStop < t0 || trStart > t1
                continue
            end

            x0 = max(trStart, t0);
            x1 = min(trStop, t1);

            if trialSpectStruct.validTrial(tr)
                switch trialSpectStruct.trialStimCode(tr)
                    case 0
                        col = [0.5 0.5 0.5];          % baseline gray
                    case 1
                        col = [0 0.4470 0.7410];      % blue
                    case 2
                        col = [0.9290 0.6940 0.1250]; % yellow
                    case 3
                        col = [0.4660 0.6740 0.1880]; % green
                    case 4
                        col = [1      0.5333 0]; % orange
                    case 5
                        col = [0.4940 0.1840 0.5560]; % purple
                    otherwise
                        col = [1 0.5 1];            % fallback magenta
                end
            elseif ~trialSpectStruct.validTrial(tr)
                col = [1 0 0];                        % invalid trial = red
            else
                continue
            end

            patch(ax, ...
                [x0 x1 x1 x0], ...
                [yl(1) yl(1) yl(2) yl(2)], ...
                col, ...
                'FaceAlpha', 0.18, ...
                'EdgeColor', 'none', ...
                'HandleVisibility', 'off', ...
                'Tag', 'TrialPatch');
        end

        uistack(findobj(ax, 'Type', 'line'), 'top')
    end

%% Trial line helper
    function plotTrialLines(ax, t0, t1)
        for tr = 1:numel(trialSpectStruct.trialStimCode)
            trStart = trialSpectStruct.trialStartIdx(tr) / fs;
            trStop  = trialSpectStruct.trialStopIdx(tr) / fs;

            if trStart >= t0 && trStart <= t1
                xline(ax, trStart, ...
                    'Color', [0 0.4470 0.7410], ...
                    'LineWidth', 2, ...
                    'HandleVisibility', 'off');
            end

            if trStop >= t0 && trStop <= t1
                xline(ax, trStop, ...
                    'r', ...
                    'LineWidth', 2, ...
                    'HandleVisibility', 'off');
            end
        end
    end

%% Copy Axes helper
    function copyAxesContent(sourceAx, targetAx)
        % Copy graphics children
        children = allchild(sourceAx);
        copyobj(flipud(children), targetAx);

        % Copy axis properties
        targetAx.XLim = sourceAx.XLim;
        targetAx.YLim = sourceAx.YLim;
        targetAx.XScale = sourceAx.XScale;
        targetAx.YScale = sourceAx.YScale;
        targetAx.CLim = sourceAx.CLim;

        % Copy labels
        xlabel(targetAx, sourceAx.XLabel.String);
        ylabel(targetAx, sourceAx.YLabel.String);

        % Copy color/axis direction
        set(targetAx, ...
            'YDir', sourceAx.YDir, ...
            'Layer', sourceAx.Layer, ...
            'Box', sourceAx.Box);
    end
end
   