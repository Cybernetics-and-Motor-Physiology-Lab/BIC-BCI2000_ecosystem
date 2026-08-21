function classifierStruct = parseBCI2000Classifier(classifierParam)
% parseBCI2000Classifier Parse BCI2000 classifier parameter Value field.
%
% INPUT
%   classifierParam : struct containing field .Value, expected as N x 4 cell
%                     Columns:
%                       1 - input channel
%                       2 - input element (bin / time / frequency)
%                       3 - output channel
%                       4 - weight
%
% OUTPUT
%   classifierStruct : struct with one field per output channel
%                     Each field contains:
%                       .inputChannelRaw      - original channel entries
%                       .inputChannelNumeric  - numeric channel indices where possible
%                       .inputElementRaw      - original bin entries
%                       .inputElementNumeric  - numeric bin values where possible
%                       .weights              - numeric weights
%                       .binSpan              - [min max] over numeric bins
%
% EXAMPLE
%   S = parseBCI2000Classifier(parameters.Classifier);
%   S.Out1.binSpan
%
% NOTES
%   - Output field names are converted into valid MATLAB field names.
%   - Non-numeric input elements/channels that cannot be parsed are returned
%     as NaN in the corresponding numeric arrays.

    % Validate input
    if ~isstruct(classifierParam) || ~isfield(classifierParam, 'Value')
        error('Input must be a struct containing field "Value".');
    end

    val = classifierParam.Value;

    if ~iscell(val) || size(val,2) ~= 4
        error('classifierParam.Value must be an N x 4 cell array.');
    end

    nRows = size(val,1);

    inputChRaw   = val(:,1);
    inputElemRaw = val(:,2);
    outputChRaw  = val(:,3);
    weightRaw    = val(:,4);

    % Normalize output channel labels for grouping
    outputLabels = cell(nRows,1);
    for i = 1:nRows
        outputLabels{i} = normalizeLabel(outputChRaw{i}, 'Out');
    end

    uniqueOutputs = unique(outputLabels, 'stable');
    classifierStruct = struct();

    for k = 1:numel(uniqueOutputs)
        thisOut = uniqueOutputs{k};
        idx = strcmp(outputLabels, thisOut);

        inChRaw_k   = inputChRaw(idx);
        inElemRaw_k = inputElemRaw(idx);
        wRaw_k      = weightRaw(idx);

        inChNum   = nan(sum(idx),1);
        inElemNum = nan(sum(idx),1);
        weights   = nan(sum(idx),1);

        for j = 1:sum(idx)
            inChNum(j)   = parseNumericToken(inChRaw_k{j});
            inElemNum(j) = parseNumericToken(inElemRaw_k{j});
            weights(j)   = parseNumericToken(wRaw_k{j});
        end

        fieldName = matlab.lang.makeValidName(thisOut);

        classifierStruct.(fieldName).inputChannelRaw = inChRaw_k;
        classifierStruct.(fieldName).inputChannelNumeric = inChNum;
        classifierStruct.(fieldName).inputElementRaw = inElemRaw_k;
        classifierStruct.(fieldName).inputElementNumeric = inElemNum;
        classifierStruct.(fieldName).weights = weights;

        validBins = inElemNum(~isnan(inElemNum));
        if isempty(validBins)
            classifierStruct.(fieldName).binSpan = [NaN NaN];
        else
            classifierStruct.(fieldName).binSpan = [min(validBins) max(validBins)];
        end
    end
end


function label = normalizeLabel(x, prefix)
% Convert output-channel entry into a stable string label.

    if nargin < 2
        prefix = 'Label';
    end

    if isnumeric(x)
        label = sprintf('%s%d', prefix, x);
        return;
    end

    if isstring(x)
        x = char(x);
    end

    if ischar(x)
        x = strtrim(x);
        numVal = str2double(x);
        if ~isnan(numVal)
            label = sprintf('%s%d', prefix, round(numVal));
        else
            label = x;
        end
        return;
    end

    label = sprintf('%sUnknown', prefix);
end


function val = parseNumericToken(x)
% Extract numeric value from mixed input like:
%   3, '3', 'Ch3', '120ms', '120 ms', '15.5', etc.

    if isnumeric(x) && isscalar(x)
        val = double(x);
        return;
    end

    if isstring(x)
        x = char(x);
    end

    if ischar(x)
        tok = regexp(strtrim(x), '[-+]?\d*\.?\d+', 'match', 'once');
        if isempty(tok)
            val = NaN;
        else
            val = str2double(tok);
        end
        return;
    end

    val = NaN;
end