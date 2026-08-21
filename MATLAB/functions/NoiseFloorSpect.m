function [mps] = NoiseFloorSpect(data,srate,freqs)    
    % calculate power spectra 
    disp('calculating power spectra')
    num_chans=size(data,2);
    win_fxn=hann(2*srate); % windowing function
    noverlap=floor(length(win_fxn) * 0.5); %overlap between calculations
    
    mps = zeros(num_chans,length(freqs));
    for chan=1:num_chans    
        [mps(chan,:),freqs] = pwelch(data(:,chan),win_fxn,noverlap,freqs,srate); % mean power spectrum across whole experiment
%         [mps(chan,:),f] = plomb(data(:,chan),srate,f,'psd'); % mean power spectrum across whole experiment
%         [ps(chan,:,curr_trial),f] = plomb(data(find(trialnr==curr_trial),chan),srate,f,'psd'); % single trial power spectra
    end
    
    % normalize power spectral densities
    % ps_mean=mean(ps(:,ccep_p>0),2);
  end