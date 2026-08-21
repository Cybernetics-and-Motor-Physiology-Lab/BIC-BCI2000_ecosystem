function T = parseBCI2000_parameters_values(parameters, paramField)
% parseBCI2000_parameters_values Convert BCI2000 parameter field into usable MATLAB table.
% DESCRIPTION
%   This function converts a BCI2000 parameter field (as returned by
%   load_bcidat) into a MATLAB table with numerically usable values.
%
%   BCI2000 parameters are stored as structured objects containing multiple
%   representations of the same data:
%   - Value         : cell array (human-readable, possibly nested)
%   - NumericValue  : numeric array (machine-readable, may contain NaNs)
%   - Type          : parameter type (e.g., 'int', 'float', 'matrix')
%   - RowLabels     : optional row names
%   - ColumnLabels  : optional column names
%
%   This function merges Value and NumericValue into a consistent table:
%   - Numeric values are taken from NumericValue when valid
%   - Missing numeric entries (NaNs) are reconstructed from Value
%   - Nested cells are converted into numeric vectors
%   - Expressions (e.g., 'StimulusCode==1') are preserved as strings
%
% INPUTS
%   parameters : struct
%       Structure returned by load_bcidat, containing parameter fields.
%   paramField : char or string
%       Name of the parameter field to extract (e.g., 'StimTriggers').
%
% BCI2000 PARAMETER STRUCTURE
%   parameters.(paramField) contains:
%   .Type           : string describing type ('int', 'float', 'matrix', etc.)
%   .Value          : cell array with following structure:
%       - rows 1:end = data rows
%       Each cell entry can be:
%           1) Numeric string:
%               '250', '10', '-3.5'
%           2) Numeric string with units:
%               '250ms', '10mA', '5Hz'
%           3) Nested cell array:
%               {'25','26','27'}
%           4) Expression / logical condition:
%               'StimulusCode==1', 'Feedback==0'
%   .NumericValue   : numeric array (same shape as Value data)
%       - contains parsed numeric values
%       - contains NaN where parsing failed
%   .RowLabels      : optional cell array of row names
%   .ColumnLabels   : optional cell array of column names
%
% OUTPUT
%   T : MATLAB table
%       Each column corresponds to a parameter variable.
%       Column types:
%           - Numeric columns:
%               when all values are scalar numeric
%           - Cell columns
%               when values contain vectors, strings, or mixed types
%       Nested values become numeric vectors stored in cells.
%       Expressions remain as string entries.
%
% CONVERSION RULES
%   For each cell entry:
%   1) If NumericValue is valid (not NaN):
%        → use NumericValue
%   2) If Value contains numeric string:
%        '250' → 250
%        '250ms' → 250
%   3) If Value contains nested cell:
%        {'1','2','3'} → [1 2 3]
%   4) If Value contains logical/parameter expression:
%        'StimulusCode==1'
%        → kept as string
%   5) If entry cannot be interpreted:
%        → warning issued
%        → value set to NaN
% EXAMPLE
%   T = parseBCI2000_parameters_values(parameters, 'StimTriggers');
%   Access:
%       T.Amplitude        % numeric column
%       T.PulseTimes{1}    % numeric vector from nested cell
%       T.Condition        % string condition
% AUTHOR
%   Frederik Lampert, Mayo Clinic 2026

    % --- Input validation ---
    if ~isstruct(parameters) || isempty(parameters)
        error('parameters must be a nonempty struct.');
    end

    if ~(ischar(paramField) || (isstring(paramField) && isscalar(paramField)))
        error('paramField must be a char vector or string scalar.');
    end

    paramField = char(paramField);
    if ~isfield(parameters, paramField)
        error('Field "%s" not found in parameters.', paramField);
    end

    P = parameters.(paramField);
    requiredFields = {'Type','Value','NumericValue'};
    for i = 1:numel(requiredFields)
        if ~isfield(P, requiredFields{i})
            error('parameters.%s must contain field "%s".', paramField, requiredFields{i});
        end
    end

    valueCell = P.Value;
    numericValue = P.NumericValue;
    if ~iscell(valueCell) || isempty(valueCell)
        error('parameters.%s.Value must be a nonempty cell array.', paramField);
    end

    % --- Column names ---
    if isfield(P, 'ColumnLabels') && ~isempty(P.ColumnLabels)
        varNames = normalizeCellstr(P.ColumnLabels);
        varNames = cellfun(@localToVarName, varNames, 'UniformOutput', false);
    else
        % Assign generic names
        nCols = size(valueCell, 2);
        varNames = arrayfun(@(i) sprintf('Var%d', i), 1:nCols, 'UniformOutput', false);
    end
    dataValueCell = valueCell;
    nRows = size(dataValueCell, 1);
    nCols = size(dataValueCell, 2);

    % --- Row names ---
    useRowNames = false;
    rowNames = {};
    if isfield(P, 'RowLabels') && ~isempty(P.RowLabels)
        rowNames = normalizeCellstr(P.RowLabels);
        if numel(rowNames) == nRows
            rowNames = matlab.lang.makeUniqueStrings(rowNames);
            useRowNames = true;
        end
    end

    % --- Numeric alignment ---
    numericCell = cell(nRows, nCols);
    if isnumeric(numericValue) && isequal(size(numericValue), [nRows nCols])
        for r = 1:nRows
            for c = 1:nCols
                numericCell{r,c} = numericValue(r,c);
            end
        end
    else
        numericCell(:) = {NaN};
    end

    % --- Parse entries ---
    parsed = cell(nRows, nCols);
    for r = 1:nRows
        for c = 1:nCols
            nv = numericCell{r,c};
            if isnumeric(nv) && isscalar(nv) && ~isnan(nv)
                parsed{r,c} = nv;
            else
                parsed{r,c} = parseValueEntry(dataValueCell{r,c}, r, c);
            end
        end
    end

    % --- Convert to table ---
    tableCols = cell(1, nCols);
    for c = 1:nCols
        col = parsed(:,c);
        if all(cellfun(@(x) isnumeric(x) && isscalar(x), col))
            tableCols{c} = cell2mat(col);
        else
            tableCols{c} = col;
        end
    end

    if useRowNames
        T = table(tableCols{:}, 'VariableNames', varNames, 'RowNames', rowNames);
    else
        T = table(tableCols{:}, 'VariableNames', varNames);
    end
end

%% Helper functions
function out = parseValueEntry(x, r, c)
    if isempty(x)
        out = NaN;
        return
    end

    if isnumeric(x)
        out = x;
        return
    end

    if ischar(x) || (isstring(x) && isscalar(x))
        s = strtrim(char(string(x)));
        % Try numeric extraction first, including units like '250ms' or '10mA'
        numVal = parseNumericString(s);
        if ~isnan(numVal)
            out = numVal;
        else
            % Keep nonnumeric expressions, e.g. 'StimulusCode==1'
            out = string(s);
        end
        return
    end

    if iscell(x)
        vals = flattenNumericCell(x);
        if isempty(vals)
            warning('Could not convert nested cell at row %d, column %d.', r, c);
            out = NaN;
        else
            out = vals;
        end
        return
    end

    warning('Could not convert entry at row %d, column %d. Leaving as NaN.', r, c);
    out = NaN;
end

function val = parseNumericString(s)
% Extract numeric value from strings like:
% '250', '250ms', '10mA', '-2.5', '1e-3'.
% Returns NaN for expressions like 'StimulusCode==1'.
    if contains(s, {'==','~=','>=','<=','>','<','&&','||'})
        val = NaN;
        return
    end

    tok = regexp(s, '^\s*([-+]?(?:\d+\.?\d*|\.\d+)(?:[eE][-+]?\d+)?)', ...
        'tokens', 'once');
    if isempty(tok)
        val = NaN;
    else
        val = str2double(tok{1});
    end
end

function vals = flattenNumericCell(x)
% Recursively flatten nested numeric/string-numeric cells into a row vector.
    vals = [];
    for k = 1:numel(x)
        thisVal = x{k};
        if isempty(thisVal)
            continue
        elseif isnumeric(thisVal)
            vals = [vals, double(thisVal(:).')]; %#ok<AGROW>
        elseif ischar(thisVal) || (isstring(thisVal) && isscalar(thisVal))
            numVal = parseNumericString(char(string(thisVal)));
            if ~isnan(numVal)
                vals = [vals, numVal]; %#ok<AGROW>
            end
        elseif iscell(thisVal)
            vals = [vals, flattenNumericCell(thisVal)]; %#ok<AGROW>
        end
    end
end

function names = normalizeCellstr(x)
    if isstring(x)
        names = cellstr(x);
    elseif ischar(x)
        names = cellstr(string(x));
    elseif iscell(x)
        names = cellfun(@(v) char(string(v)), x, 'UniformOutput', false);
    else
        error('Labels must be char, string, or cell array.');
    end
    names = reshape(names, 1, []);
end

function name = localToVarName(x)
    if isstring(x)
        x = char(x);
    end

    if ~ischar(x)
        x = char(string(x));
    end
    name = matlab.lang.makeValidName(strtrim(x));
end