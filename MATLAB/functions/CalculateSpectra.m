function [f,data_epoch_spectra,data_epoch] = ...
    CalculateSpectra(data_epoch,wndw,epch_t,ovrlp,fs,frqs)

% Function to calculate powerspectrum for across several ECoG channels and
% epochs using Matlab's pwelch function.
% 
% Input:
% data_epoch: channels X epochs X time
% wndw - window size
% fft_t - time window of epoch
%
%   
%
% Output:
% data_epoch_spectra is channels X epochs X frequencies
%
% Example:
% [f,data_epoch_spectra,data_epoch] = ...
%     ecog_spectra(data_epoch,stims,fft_w,fft_t,fft_ov,srate,reg_erp)
%
%
% Created by Dora Hermes, 2017, moddified by Frederik Lampert 2024

% 
% % calculate spectra to get length of f to initialize spectra
% [~,f] = pwelch(squeeze(data_epoch(1,1:fft_t,1)),fft_w,fft_ov,srate,srate);
% data_epoch_spectra = zeros(size(data_epoch,1),size(data_epoch,2),length(f));
% clear Pxx

% calculate powerspectra
for k = 1:size(data_epoch,1)%channels
    disp(['fft el ' int2str(k)])
    x = squeeze(data_epoch(k,:,:))';
    [Pxx,f] = pwelch(x(fft_t,:),fft_w,fft_ov,srate,srate);
    data_epoch_spectra(k,:,:) = Pxx';
end