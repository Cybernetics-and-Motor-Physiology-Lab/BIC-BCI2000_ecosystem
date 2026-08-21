function [recon_sig] = interpolate_lost_samples(raw_sig, sig_states, method)
    if nargin < 3
        method = 'linear';
    end

    recon_sig = raw_sig;            % Copy the orginal signal
    recon_sig(sig_states.ImplantLostSample==1,:) = NaN; % Drop the lost values
    recon_sig = fillmissing(recon_sig,method);    % Fill the missing values by linear interpolation
end
