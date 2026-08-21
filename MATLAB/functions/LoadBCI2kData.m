function [signal, states, parameters, total_samples, file_path] = LoadBCI2kData(filePath, interpolationMethod)
% LoadBCI2kData Load calibrated BCI2000 .dat data and post-process states.
% Dependency : BCI2000 mex tools https://www.bci2000.org/mediawiki/index.php/User_Reference:Matlab_MEX_Files
% Syntax
%   [signal, states, parameters, total_samples, file_path] = LoadBCI2kData()
%   [signal, states, parameters, total_samples, file_path] = LoadBCI2kData(filePath)
%   [signal, states, parameters, total_samples, file_path] = LoadBCI2kData(filePath, interpolationMethod)
%
% Description
%   This function loads a BCI2000 .dat file using LOAD_BCIDAT with the
%   '-calibrated' option, optionally interpolates samples marked as lost by
%   the ImplantLostSample state, and converts selected integer-encoded
%   states to single precision.
%
%   If filePath is not provided or is empty, the user is prompted to select
%   a .dat file through a GUI.
%
% Inputs
%   filePath            - Optional. Character vector or string scalar
%                         specifying the full path to a .dat file.
%
%   interpolationMethod - Optional. Missing-sample interpolation method.
%                         Must be one of:
%                           'none'     - leave NaNs as-is
%                           'previous' - previous non-missing value
%                           'next'     - next non-missing value
%                           'nearest'  - nearest non-missing value
%                           'linear'   - linear interpolation
%                           'spline'   - cubic spline interpolation
%                           'pchip'    - shape-preserving cubic interpolation
%                           'makima'   - modified Akima interpolation
%                         Default: 'linear'
%
% Outputs
%   signal        - Numeric matrix [samples x channels], calibrated signal.
%   states        - Struct containing BCI2000 states. Matching SGNL*,
%                   NormalizerOffset%d, and NormalizerGain%d fields are
%                   typecast to single precision.
%   parameters    - Struct returned by LOAD_BCIDAT containing file metadata.
%   total_samples - Total number of samples returned by LOAD_BCIDAT.
%   file_path     - Full path to the loaded .dat file.
%
% Notes
%   - Lost samples are identified using states.ImplantLostSample when that
%     field exists.
%   - If interpolationMethod is 'none', lost samples are left as NaN.
%   - This function assumes selected state fields store single-precision
%     values encoded in an integer container and uses TYPECAST to recover
%     the original single values.
%
% Example
%   [signal, states, parameters, total_samples, file_path] = LoadBCI2kData();
%
%   [signal, states] = LoadBCI2kData('/path/to/file.dat', 'pchip');


    % Input validation and defaults
    if nargin < 1
        filePath = [];
    end
    if nargin < 2 || isempty(interpolationMethod)
        interpolationMethod = 'linear';
    end

    validInterp = {'none','previous','next','nearest','linear','spline','pchip','makima'};

    if ~(ischar(interpolationMethod) || (isstring(interpolationMethod) && isscalar(interpolationMethod)))
        error('interpolationMethod must be a character vector or string scalar.');
    end
    interpolationMethod = char(string(interpolationMethod));
    if ~ismember(interpolationMethod, validInterp)
        error('Invalid interpolationMethod "%s". Must be one of: %s.', ...
            interpolationMethod, strjoin(validInterp, ', '));
    end

    if ~(isempty(filePath) || ischar(filePath) || (isstring(filePath) && isscalar(filePath)))
        error('filePath must be empty, a character vector, or a string scalar.');
    end
    filePath = char(string(filePath));

    % --------------------------
    % Resolve file path
    % --------------------------
    % User may customize this path; otherwise fallback to current folder
    default_pth = ['/Users/lampert.frederik/Documents/Mayo/CorTec/data/' filesep];
    if ~isfolder(default_pth)
        default_pth = [pwd filesep];
    end

    if isempty(filePath)
        [file, folder] = uigetfile(fullfile(default_pth, '*.dat'), 'Select a File');
        if isequal(file, 0)
            error('No file was selected.');
        end
        file_path = fullfile(folder, file);
    else
        if ~isfile(filePath)
            error('The specified file does not exist: %s', filePath);
        end
        [~, ~, ext] = fileparts(filePath);
        if ~strcmpi(ext, '.dat')
            error('filePath must point to a .dat file.');
        end
        file_path = filePath;
    end

    % --------------------------
    % Load data
    % --------------------------
    [signal, states, parameters, total_samples] = load_bcidat(file_path, '-calibrated');

    % Output validation after load
    if ~isnumeric(signal) || ndims(signal) ~= 2
        error('Loaded signal must be a 2D numeric matrix.');
    end
    if ~isstruct(states)
        error('Loaded states must be a struct.');
    end
    if ~isstruct(parameters)
        error('Loaded parameters must be a struct.');
    end
    if ~isscalar(total_samples) || ~isnumeric(total_samples)
        error('Loaded total_samples must be a numeric scalar.');
    end

    % --------------------------
    % Mark and interpolate lost samples
    % --------------------------
    if isfield(states, 'ImplantLostSample')
        lostMask = states.ImplantLostSample == 1;
        if ~isvector(lostMask)
            error('states.ImplantLostSample must be a vector.');
        end
        lostMask = lostMask(:);
        if numel(lostMask) ~= size(signal, 1)
            error('states.ImplantLostSample length (%d) does not match signal samples (%d).', ...
                numel(lostMask), size(signal, 1));
        end
        signal(lostMask, :) = NaN;
        if ~strcmp(interpolationMethod, 'none')
            signal = fillmissing(signal, interpolationMethod);
        end
    end

    % --------------------------
    % Convert selected states from integer payload to single
    % --------------------------
    fn = fieldnames(states);
    % Match:
    %   1) any field starting with 'SGNL'
    %   2) NormalizerOffset%d
    %   3) NormalizerGain%d
    isSGNL = startsWith(fn, 'SGNL');
    isNorm = ~cellfun(@isempty, regexp(fn, '^Normalizer(?:Offset|Gain)\d+$', 'once'));
    fieldsToConvert = fn(isSGNL | isNorm);

    for k = 1:numel(fieldsToConvert)
        thisField = fieldsToConvert{k};
        rawVal = states.(thisField);
        if ~isnumeric(rawVal)
            warning('Skipping field "%s": expected numeric input for typecast.', thisField);
            continue;
        end

        % Typecast elementwise while preserving original size
        try
            convertedVal = typecast(rawVal(:), 'single');
            if numel(convertedVal) == numel(rawVal)
                states.(thisField) = reshape(convertedVal, size(rawVal));
            else
                warning('Skipping field "%s": typecast size mismatch.', thisField);
            end
        catch ME
            warning('Skipping field "%s": %s', thisField, ME.message);
        end
    end
end