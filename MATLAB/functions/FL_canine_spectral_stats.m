function rvals=FL_canine_spectral_stats(nps,f,tr_sc,LFB_range, HFB_range, BF_correct)
% FL_canine_spectral_stats - Compute spectral statistics for broadband and low-frequency ranges.
%
% Syntax:
%   rvals = FL_canine_spectral_stats(nps, f, tr_sc, LFB_range, HFB_range, BF_correct)
%
% Description:
%   This function calculates statistical metrics for normalized trial power spectra (nps).
%   It compares broadband frequency activity (HFB) and low-frequency oscillations (LFB)
%   between two trial conditions (tr_sc==1 vs. tr_sc==0) using signed r-squared and
%   unpaired t-tests. Bonferroni correction can optionally be applied to p-values.
%
% Inputs:
%   nps        - (NxFxT array) Normalized trial power spectra, where:
%                N = number of channels
%                F = number of frequencies
%                T = number of trials
%   f          - (vector) Frequency range corresponding to the power spectra.
%   tr_sc      - (vector) Trial stimulus code indicating trial type (e.g., 1 or 0).
%   LFB_range  - (vector) Indices specifying the low-frequency oscillations range in `f`.
%   HFB_range  - (vector) Indices specifying the broadband frequency range in `f`.
%   BF_correct - (logical) Whether to apply Bonferroni correction to p-values (true/false).
%
% Outputs:
%   rvals      - (struct) Structure containing results:
%                - HFB_range: Indices for broadband frequency range.
%                - LFB_range: Indices for low-frequency oscillations range.
%                - HFB_trials: Mean power for HFB range across trials.
%                - LFB_trials: Mean power for LFB range across trials.
%                - r_HFB: Signed r-squared values for HFB trials.
%                - r_LFB: Signed r-squared values for LFB trials.
%                - p_HFB: P-values for HFB trials (unpaired t-test).
%                - p_LFB: P-values for LFB trials (unpaired t-test).
%                - rmap: Signed r-squared feature map across frequencies.
%
% Example:
%   % Example inputs
%   nps = rand(10, 100, 50); % 10 channels, 100 frequencies, 50 trials
%   f = 1:100;               % Frequency range 1 to 100 Hz
%   tr_sc = [ones(1, 25), zeros(1, 25)]; % 25 trials for each condition
%   LFB_range = 1:10;        % Low frequencies (1-10 Hz)
%   HFB_range = 50:100;      % Broadband frequencies (50-100 Hz)
%   BF_correct = true;       % Apply Bonferroni correction
%
%   % Compute spectral statistics
%   rvals = FL_canine_spectral_stats(nps, f, tr_sc, LFB_range, HFB_range, BF_correct);
%
% Notes:
%   - Ensure `nps` has the correct dimensions (channels, frequencies, trials).
%   - The `rsa` function must be defined separately for computing signed r-squared.
%   - If Bonferroni correction is enabled, p-values are adjusted based on the number of channels.
%
% See also:
%   ttest2
%
% Author:
%   Kai J Miller, Mayo Clinic , 2024
% Modified by:
%   Frederik Lampert, Mayo Clinic, 2025
% Documentation creted using ChatGPT4o


num_chans=size(nps,1); % Number of channels

% Select average broadband and low-frequency trials spectra
HFB_trials=squeeze(mean(nps(:,HFB_range,:),2)); % Select broadband frequency range
LFB_trials=squeeze(mean(nps(:,LFB_range,:),2)); % Select low-frequency oscilations range 

% Initialize output variables
r_HFB = zeros(1, num_chans);
r_LFB = zeros(1, num_chans);
p_HFB = zeros(1, num_chans);
p_LFB = zeros(1, num_chans);
rmap = zeros(num_chans, length(f));

% Perform comparisons for each channel
for chan = 1:num_chans
    % Compute signed r-squared for HFB and LFB trials
    r_HFB(chan) = rsa(HFB_trials(chan, tr_sc == 1), HFB_trials(chan, tr_sc == 0));
    r_LFB(chan) = rsa(LFB_trials(chan, tr_sc == 1), LFB_trials(chan, tr_sc == 0));

    % Perform unpaired t-tests for HFB and LFB trials
    [~, p_HFB(chan)] = ttest2(HFB_trials(chan, tr_sc == 1), HFB_trials(chan, tr_sc == 0));
    [~, p_LFB(chan)] = ttest2(LFB_trials(chan, tr_sc == 1), LFB_trials(chan, tr_sc == 0));

    % Apply Bonferroni correction if specified
    if BF_correct
        p_HFB(chan) = p_HFB(chan) * num_chans;
        p_LFB(chan) = p_LFB(chan) * num_chans;
    end

    % Compute feature map (signed r-squared for all frequencies)
    for freq = 1:length(f)
        rmap(chan, freq) = rsa(nps(chan, freq, tr_sc == 1), nps(chan, freq, tr_sc == 0));
    end
end

% Wrap results into a structure
rvals.HFB_range = HFB_range;
rvals.LFB_range = LFB_range;
rvals.HFB_trials = HFB_trials;
rvals.LFB_trials = LFB_trials;
rvals.r_HFB = r_HFB;
rvals.r_LFB = r_LFB;
rvals.p_HFB = p_HFB;
rvals.p_LFB = p_LFB;
rvals.rmap = rmap;

end
