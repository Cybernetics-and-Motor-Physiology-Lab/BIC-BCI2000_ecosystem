function [cpsdMatrix, freqVec] = channelCrossPSD(signal,fs, window, overlap, nfft)
% channelCrossPSD Compute the cross power spectral density (CPSD) matrix between channels.
%
%   cpsdMatrix = channelCrossPSD(signal)
%   cpsdMatrix = channelCrossPSD(signal, fs, window, overlap, nfft)
%
% INPUTS:
%   signal  - Input signal matrix [timeSamples x channels].
%   fs      - (Optional) Sampling frequency in Hz. Default is 1000 Hz.
%   window  - (Optional) Analysis window (default: periodic Hann window of 2048 samples).
%   overlap - (Optional) Number of overlapping samples between windows (default: 50% overlap).
%   nfft    - (Optional) Number of FFT points (default: length of window).
%
% OUTPUT:
%   cpsdMatrix - 3D complex array [channels x channels x frequencies], where each slice
%                (cpsdMatrix(:,:,f)) contains the cross power spectral density between all
%                channel pairs at frequency bin f.
%   freqVec - Frequencies at which the STFT is evaluated, returned as a vector.
%
% DESCRIPTION:
%   This function computes the Cross Power Spectral Density (CPSD) between all pairs
%   of channels in the input signal using the Short-Time Fourier Transform (STFT).
%   It assumes the signal is real-valued and returns a Hermitian matrix at each frequency.
%
% NOTES:
%   - If the input signal is [channels x timeSamples], it will be automatically transposed.
%   - Only the upper triangle of the matrix is directly computed; the lower triangle
%     is filled using Hermitian symmetry.
%   - Normalization accounts for sampling frequency and window energy.
%
% EXAMPLE USAGE:
%   [cpsdMat, freqs] = channelCrossPSD(mySignal, 2000, hann(1024), 512, 1024);
%
% AUTHOR:
%   Frederik Lampert, Mayo Clinic 2025

% Default arguments
if nargin < 5, nfft = length(window); end
if nargin < 4, overlap = floor(length(window)/2); end
if nargin < 3, window = hann(2048, 'periodic'); end
if nargin < 2, fs = 1000; end

% Ensure [timeSamples x channels]
if size(signal, 2) > size(signal,1)
    signal = signal';
end

% STFT computation
[STFTs,freqVec] =  stft(signal,fs, ...
    'Window', window, ...
    'OverlapLength', overlap, ...
    'FFTLength', nfft, ...
    'FrequencyRange', 'onesided');

% Precompute constants
nCh = size(signal, 2); 
cpsdMatrix = complex(nan(nCh, nCh, length(freqVec)));  % Preallocate complex Cross-PSD matrix
winNorm = fs * sum(window.^2); 

% CPSD computation
for ch1 = 1:nCh
    for ch2 = ch1:nCh  % Only upper triangle (ch1 <= ch2)
        % Cross-PSD at each time-frequency bin
        cpsdMatrix(ch1, ch2, :) = mean(STFTs(:,:,ch1) .* conj(STFTs(:,:,ch2)) / winNorm,2);        
        if ch1 ~= ch2
            cpsdMatrix(ch2, ch1, :) = conj(cpsdMatrix(ch1, ch2, :));  % Hermitian symmetry
        end
    end
end


%% Old (unefficent) code
% cpsdMatrix = nan(nCh, nCh, nFreq);  % Preallocate
% parfor ch1 = 1:nCh
%     localMatrix = nan(nCh, nFreq);  % Each worker writes locally
%     sig1 = signal(:, ch1);
%     for ch2 = ch1+1:nCh
%         sig2 = signal(:, ch2);
%         [Cxy, ~] = cpsd(sig1, sig2, window, overlap, freqVec, fs);
%         localMatrix(ch2, :) = Cxy;
%     end
%     cpsdMatrix(ch1, :, :) = localMatrix;
% end
% 
% % Fill lower triangle (because parfor can't do it inside easily)
% for ch1 = 1:nCh
%     for ch2 = ch1+1:nCh
%         cpsdMatrix(ch2, ch1, :) = cpsdMatrix(ch1, ch2, :);
%     end
% end