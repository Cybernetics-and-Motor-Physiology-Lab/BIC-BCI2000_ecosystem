function [outputSignal, outputLabels, info] = applyBCI2kSpatialFilter(signal, parameters, channelLabels)
% applyBCI2kSpatialFilter Apply a BCI2000 SpatialFilter to signal data.
%
% Syntax
%   [outputSignal, outputLabels, info] = applyBCI2kSpatialFilter(signal, parameters)
%   [outputSignal, outputLabels, info] = applyBCI2kSpatialFilter(signal, parameters, channelLabels)
%
% Inputs
%   signal        - Numeric matrix [samples x channels].
%   parameters    - BCI2000 parameters struct. Must contain at least:
%                     parameters.SpatialFilterType.NumericValue
%                   and, depending on type:
%                     parameters.SpatialFilter
%                     parameters.SpatialFilterCAROutput
%                     parameters.SpatialFilterMissingChannels
%   channelLabels - Optional string/cellstr labels for input channels.
%                   Default: "Ch1", "Ch2", ..., "ChN"
%
% Outputs
%   outputSignal  - Filtered signal [samples x outputChannels].
%   outputLabels  - Output channel labels as string row vector.
%   info          - Struct with metadata about the applied filter.
%
% Filter types
%   0: none
%   1: full matrix
%   2: sparse matrix
%   3: common average reference (CAR)
%
% Notes
%   - Full matrix: output = signal * W', where W is [outCh x inCh].
%   - Sparse matrix: SpatialFilter rows are interpreted as [input, output, weight].
%   - CAR: average is computed over all input channels and subtracted from
%     selected output channels or all channels when SpatialFilterCAROutput is empty.

    %% Validate inputs
    narginchk(2, 3);

    validateattributes(signal, {'numeric'}, {'2d','nonempty','real'}, mfilename, 'signal', 1);
    if size(signal,2) > size(signal,1)
        % keep the user's intended convention explicit
        warning('Input signal has more columns than rows. Assuming [samples x channels] was intended and leaving unchanged.');
    end

    if ~isstruct(parameters)
        error('parameters must be a struct.');
    end

    nInCh = size(signal, 2);

    if nargin < 3 || isempty(channelLabels)
        channelLabels = "CH" + (0:nInCh-1);
    else
        if iscell(channelLabels)
            channelLabels = string(channelLabels);
        end
        if ~(isstring(channelLabels) || ischar(channelLabels))
            error('channelLabels must be a string array, char array, or cellstr.');
        end
        channelLabels = string(channelLabels);
        if numel(channelLabels) ~= nInCh
            error('Number of channelLabels (%d) must match number of input channels (%d).', numel(channelLabels), nInCh);
        end
    end
    channelLabels = reshape(channelLabels, 1, []);

    %% Extract required parameters
    mustHaveField(parameters, 'SpatialFilterType');
    filterType = extractNumericScalar(parameters.SpatialFilterType, 'SpatialFilterType');

    % Optional / type-dependent
    spatialFilterRaw = [];
    if isfield(parameters, 'SpatialFilter')
        spatialFilterRaw = getParamPayload(parameters.SpatialFilter);
    end

    carOutputRaw = [];
    if isfield(parameters, 'SpatialFilterCAROutput')
        carOutputRaw = getParamPayload(parameters.SpatialFilterCAROutput);
    end

    missingMode = 1; % default: Report Error, per docs
    if isfield(parameters, 'SpatialFilterMissingChannels')
        missingMode = extractNumericScalar(parameters.SpatialFilterMissingChannels, 'SpatialFilterMissingChannels');
    end

    %% Dispatch
    info = struct();
    info.FilterTypeCode = filterType;
    info.FilterTypeName = filterTypeName(filterType);
    info.InputLabels = channelLabels;

    switch filterType
        case 0  % none
            outputSignal = signal;
            outputLabels = channelLabels;

        case 1  % full matrix
            mustHaveField(parameters, 'SpatialFilter');
            W = parseFullMatrix(spatialFilterRaw);
            if size(W,2) ~= nInCh
                error('Full SpatialFilter matrix has %d input columns, but signal has %d channels.', size(W,2), nInCh);
            end

            % With MATLAB [samples x channels]: Y = X * W'
            outputSignal = signal * W.';
            outputLabels = "Out" + (1:size(W,1));

            info.SpatialFilterMatrix = W;

        case 2  % sparse matrix
            mustHaveField(parameters, 'SpatialFilter');
            sparseRows = parseSparseMatrix(spatialFilterRaw);

            [outputSignal, outputLabels, sparseMap] = applySparse(signal, sparseRows, channelLabels, missingMode);
            info.SparseRows = sparseRows;
            info.SparseMap = sparseMap;

        case 3  % CAR
            [outputSignal, outputLabels, selectedIdx] = applyCAR(signal, channelLabels, carOutputRaw, missingMode);
            info.CARSelectedChannels = selectedIdx;

        otherwise
            error('Unsupported SpatialFilterType code: %g', filterType);
    end

    info.OutputLabels = outputLabels;
end


%% ========================= Helpers =========================

function mustHaveField(S, fname)
    if ~isfield(S, fname)
        error('parameters.%s is required but missing.', fname);
    end
end

function x = extractNumericScalar(paramStruct, paramName)
    x = getParamPayload(paramStruct);
    if iscell(x)
        x = str2double(cell2mat(x));
    end
    if isempty(x) || ~isnumeric(x) || numel(x) ~= 1
        error('%s must resolve to a numeric scalar.', paramName);
    end
    x = double(x);
end

function payload = getParamPayload(p)
% Try common BCI2000 parameter struct layouts.

    if ~isstruct(p)
        payload = p;
        return;
    end

    candidateFields = {'Value','NumericValue','Matrix','Data'};
    for k = 1:numel(candidateFields)
        if isfield(p, candidateFields{k})
            payload = p.(candidateFields{k});
            return;
        end
    end

    % Fallback: if struct only has one field, use it
    fn = fieldnames(p);
    if numel(fn) == 1
        payload = p.(fn{1});
        return;
    end

    error('Could not determine payload from parameter struct.');
end

function name = filterTypeName(code)
    switch code
        case 0, name = "none";
        case 1, name = "fullMatrix";
        case 2, name = "sparseMatrix";
        case 3, name = "CAR";
        otherwise, name = "unknown";
    end
end

function W = parseFullMatrix(raw)
% Return numeric [outCh x inCh] matrix.

    if isnumeric(raw)
        W = double(raw);
        return;
    end

    if iscell(raw)
        try
            W = double(cellfun(@str2double, raw));
            return;
        catch
        end
        try
            W = cell2mat(raw);
            W = double(W);
            return;
        catch
        end
    end

    if isstring(raw) || ischar(raw)
        error('String-only full-matrix SpatialFilter payload is not supported in this parser. Provide the loaded parameter payload example if needed.');
    end

    error('Unable to parse full SpatialFilter matrix.');
end

function rows = parseSparseMatrix(raw)
% Return table-like numeric rows: [inputSpec outputSpec weight]
% First two columns may be numeric or labels in the raw BCI2000 matrix, so
% we preserve them as cell array rows.

    if isnumeric(raw)
        if size(raw,2) ~= 3
            error('Sparse SpatialFilter must have exactly 3 columns: input, output, weight.');
        end
        rows = num2cell(raw);
        return;
    end

    if iscell(raw)
        if size(raw,2) ~= 3
            error('Sparse SpatialFilter must have exactly 3 columns: input, output, weight.');
        end
        rows = raw;
        return;
    end

    error('Unable to parse sparse SpatialFilter matrix.');
end

function [Y, outLabels, mapInfo] = applySparse(X, sparseRows, inLabels, missingMode)
% sparseRows: N x 3 cell array -> {inputSpec, outputSpec, weight}

    nSamples = size(X,1);
    nRows = size(sparseRows,1);

    outSpecList = strings(1, nRows);
    validRow = false(nRows,1);
    fromIdx = nan(nRows,1);
    toIdxTmp = strings(nRows,1);
    wt = nan(nRows,1);

    % Resolve input rows and collect output labels/specs
    for r = 1:nRows
        inSpec = sparseRows{r,1};
        outSpec = sparseRows{r,2};
        weight = sparseRows{r,3};

        wt(r) = str2double(weight);

        [inIdx, inFound] = resolveChannelSpec(inSpec, inLabels);
        fromIdx(r) = inIdx;

        toIdxTmp(r) = normalizeOutputSpec(outSpec);
        outSpecList(r) = toIdxTmp(r);

        if inFound && abs(wt(r)) > eps(abs(wt(r)))
            validRow(r) = true;
        elseif ~inFound && missingMode ~= 0
            error('Sparse SpatialFilter references missing input channel "%s".', string(inSpec));
        end
    end

    % Keep only valid rows
    sparseRows = sparseRows(validRow,:);
    fromIdx = fromIdx(validRow);
    wt = wt(validRow);
    toIdxTmp = toIdxTmp(validRow);

    if isempty(sparseRows)
        Y = zeros(nSamples, 0);
        outLabels = strings(1,0);
        mapInfo = struct('fromIdx',[],'toIdx',[],'weights',[]);
        return;
    end

    % Unique output channels preserve first appearance order
    [outLabels, ~, toIdx] = unique(toIdxTmp, 'stable');
    nOut = numel(outLabels);

    Y = zeros(nSamples, nOut, 'like', X);
    for r = 1:numel(fromIdx)
        Y(:, toIdx(r)) = Y(:, toIdx(r)) + X(:, fromIdx(r)) * wt(r);
    end

    mapInfo = struct('fromIdx', fromIdx, 'toIdx', toIdx, 'weights', wt);
end

function [Y, outLabels, selectedIdx] = applyCAR(X, inLabels, carOutputRaw, missingMode)
    nSamples = size(X,1);
    nIn = size(X,2);

    avg = mean(X, 2);  % all input channels participate in CAR, per docs

    selectedIdx = [];
    outLabels = strings(1,0);

    if isempty(carOutputRaw)
        selectedIdx = 1:nIn;
        outLabels = inLabels;
    else
        specs = normalizeListPayload(carOutputRaw);
        for k = 1:numel(specs)
            [idx, found] = resolveChannelSpec(specs{k}, inLabels);
            if found
                selectedIdx(end+1) = idx; %#ok<AGROW>
                outLabels(end+1) = inLabels(idx); %#ok<AGROW>
            elseif missingMode ~= 0
                error('Invalid channel specification "%s" in SpatialFilterCAROutput(%d).', string(specs{k}), k);
            end
        end
    end

    Y = X(:, selectedIdx) - avg;
end

function [idx, found] = resolveChannelSpec(spec, inLabels)
% spec may be numeric channel index (1-based) or label string.

    found = false;
    idx = -1;

    if isnumeric(spec)
        idx = floor(double(spec));
        found = idx >= 1 && idx <= numel(inLabels);
        if ~found
            idx = -1;
        end
        return;
    end

    s = string(spec);
    s = strip(s);

    % Try numeric string first
    asNum = str2double(s);
    if ~isnan(asNum)
        idx = floor(asNum);
        found = idx >= 1 && idx <= numel(inLabels);
        if ~found
            idx = -1;
        end
        return;
    end

    % Label match
    m = find(strcmpi(inLabels, s), 1);
    if ~isempty(m)
        idx = m;
        found = true;
    end
end

function s = normalizeOutputSpec(outSpec)
% Create a stable string label for sparse-matrix output channel specs.
    if isnumeric(outSpec)
        s = "Out" + string(floor(double(outSpec)));
    else
        s = string(outSpec);
        s = strip(s);
        asNum = str2double(s);
        if ~isnan(asNum)
            s = "Out" + string(floor(asNum));
        end
    end
end

function specs = normalizeListPayload(raw)
% Normalize SpatialFilterCAROutput payload into a cell array of scalar specs.

    if isempty(raw)
        specs = {};
        return;
    end

    if isnumeric(raw)
        specs = num2cell(raw(:).');
        return;
    end

    if isstring(raw)
        specs = cellstr(raw(:)).';
        return;
    end

    if ischar(raw)
        % single token or whitespace-separated list
        parts = strsplit(strtrim(raw));
        specs = parts;
        return;
    end

    if iscell(raw)
        specs = raw(:).';
        return;
    end

    error('Unsupported SpatialFilterCAROutput payload type.');
end