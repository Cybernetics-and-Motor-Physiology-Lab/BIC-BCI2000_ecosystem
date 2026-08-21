function cohMatrix = cpsdToCoherence(cpsdMatrix)
% cpsdToCoherence Compute coherence matrix from CPSD matrix
%
% INPUT:
%   cpsdMatrix - complex CPSD matrix [channels x channels x frequencies]
%
% OUTPUT:
%   cohMatrix  - coherence matrix [channels x channels x frequencies]
%
% DESCRIPTION:
%   Computes magnitude-squared coherence:
%       Cxy(f) = |Pxy(f)|^2 / (Pxx(f)*Pyy(f))

    nCh = size(cpsdMatrix, 1);
    nFreq = size(cpsdMatrix, 3);

    cohMatrix = nan(nCh, nCh, nFreq);

    % Extract auto-spectra from diagonal
    autoPSD = zeros(nCh, nFreq);
    for ch = 1:nCh
        autoPSD(ch, :) = real(squeeze(cpsdMatrix(ch, ch, :)));
    end

    % Compute coherence
    for ch1 = 1:nCh
        for ch2 = 1:nCh
            denom = autoPSD(ch1, :) .* autoPSD(ch2, :);
            valid = denom > 0;
            coh = nan(1, nFreq);
            % Extract full CPSD vector first, then apply logical indexing
            cpsdVec = squeeze(cpsdMatrix(ch1, ch2, :)).';
            coh(valid) = abs(cpsdVec(valid)).^2 ./ denom(valid);
            cohMatrix(ch1, ch2, :) = coh;
        end
    end
end