function filtSig = fastButterFilt(signal, fs, type, stopFreq, order, verbose)
% fastButterFilt
% Apply fast Butterworth filtering using a mex backend.
%
% REQUIRED:
%   signal    : [samples x channels] or [channels x samples]
%   fs        : sampling frequency (Hz)
%   type      : 'high' | 'low' | 'notch' | 'bandpass'
%   stopFreq  : cutoff frequency (Hz)
%               - scalar for high / low / stop
%               - [f1 f2] for bandpass (f1 < f2)
%
% OPTIONAL:
%   order     : filter order (default = 3)
%   verbose   : logical, print info (default = false)

%% -------------------- Defaults --------------------
if nargin < 5 || isempty(order)
    order = 3;
end

if nargin < 6 || isempty(verbose)
    verbose = false;
end

%% -------------------- Constants --------------------
stopWidth = 1;  % Hz, used only for notch (stop) filter

%% -------------------- Input checks --------------------
validateattributes(signal, {'numeric'}, {'2d'}, mfilename, 'signal');
validateattributes(fs, {'numeric'}, {'scalar','positive'}, mfilename, 'fs');
validateattributes(order, {'numeric'}, {'scalar','integer','positive'}, mfilename, 'order');
validateattributes(verbose, {'logical','numeric'}, {'scalar'}, mfilename, 'verbose');

type = lower(type);
validTypes = {'high','low','notch','bandpass'};
assert(ismember(type, validTypes), ...
    'fastButterFilt:BadType', ...
    'Filter type must be one of: %s', strjoin(validTypes, ', '));

if strcmp(type, 'bandpass')
    assert(numel(stopFreq) == 2 && stopFreq(1) < stopFreq(2), ...
        'fastButterFilt:BadBandpassFreq', ...
        'Bandpass stopFreq must be [f1 f2] with f1 < f2');
else
    assert(isscalar(stopFreq), ...
        'fastButterFilt:BadStopFreq', ...
        'stopFreq must be scalar for %s filter', type);
end

%% -------------------- Ensure samples x channels --------------------
if size(signal,1) < size(signal,2)
    signal = signal.';
end

%% -------------------- Add mex path if needed --------------------
mexPath = '/Users/lampert.frederik/Documents/MATLAB/CnelViewer_devfork/src/mex';
if ~any(strcmp(mexPath, path))
    addpath(mexPath);
end

%% -------------------- Logging --------------------
if verbose
    fprintf('Applying %s Butterworth filter at %s Hz (order = %d)\n', ...
        type, mat2str(stopFreq), order);
end

%% -------------------- Handle NaNs --------------------
nanMask = isnan(signal);
if any(nanMask(:))
    signal(nanMask) = 0;
end

%% -------------------- Padding for edge effects --------------------
ext = min(2*fs, floor(size(signal,1)/2));

%% -------------------- Design filter --------------------
nyq = fs / 2;

switch type
    case 'high'
        Wn = stopFreq / nyq;
        [b, a] = butter(order, Wn, 'high');

    case 'low'
        Wn = stopFreq / nyq;
        [b, a] = butter(order, Wn, 'low');

    case 'notch'   % notch
        Wn = [stopFreq-stopWidth, stopFreq+stopWidth] / nyq;
        [b, a] = butter(order, Wn, 'stop');

    case 'bandpass'
        Wn = stopFreq / nyq;
        [b, a] = butter(order, Wn, 'bandpass');
end

%% -------------------- Expand coefficients per channel --------------------
nChan = size(signal,2);

bCell = cell(1, nChan);
aCell = cell(1, nChan);

for ch = 1:nChan
    % Each channel gets a CELL ARRAY of filter sections
    bCell{ch} = {b};
    aCell{ch} = {a};
end

%% -------------------- Signal extension --------------------
signalExt = wextend('ar','sym', signal, ext, 'b');

%% -------------------- Call mex filter --------------------
if isunix
    filtSig = UnixMultiThreadedFilter(bCell, aCell, signalExt);
elseif ispc
    filtSig = WinMultiThreadedFilter(bCell, aCell, signalExt);
else
    error('fastButterFilt:UnsupportedOS', ...
        'No supported mex filter for this OS.');
end

%% -------------------- Remove padding --------------------
filtSig = filtSig(ext+1:end-ext, :);

%% -------------------- Restore NaNs --------------------
if any(nanMask(:))
    filtSig(nanMask) = NaN;
end

end