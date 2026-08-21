function result = matchTimestamps(A, B, tolerance)
% matchTimestamps - Matches timestamps from vector A to vector B based on expected delay
%
% Inputs:
%   A         - Vector of timestamps when commands were sent
%   B         - Vector of timestamps when commands were received
%   tolerance - [minDelay maxDelay] in samples (e.g., [100 250])
%
% Output (struct):
%   result.pairedIdxA   - Index in B matched to each A (NaN if no match)
%   result.unmatchedA   - Indices in A with no match
%   result.multiMatchedB - Indices in B that were matched more than once
%   result.pairTable    - Table of matched A and B values (NaN if unmatched)
% Created by ChatGPT, veryfied by FL
% Nov, 2025

    if nargin < 3
        tolerance = [1 120];  % Default delay window
    end

    nA = length(A);
    nB = length(B);

    pairedIdxA = nan(nA, 1);  % For each A(i), which B(j) is matched
    matchCountB = zeros(nB, 1);  % Count how many times each B is matched

    for i = 1:nA
        diffs = B - A(i);
        inRange = diffs >= tolerance(1) & diffs <= tolerance(2);
        candidates = find(inRange);

        if ~isempty(candidates)
            [~, bestRelIdx] = min(diffs(candidates));
            bestBIdx = candidates(bestRelIdx);
            pairedIdxA(i) = bestBIdx;
            matchCountB(bestBIdx) = matchCountB(bestBIdx) + 1;
        end
    end

    unmatchedA    = find(isnan(pairedIdxA));
    unmatchedB    = find(matchCountB == 0);
    multiMatchedB = find(matchCountB > 1);

    % Optional: Table of paired timestamps
    pairTable = table(A(:), nan(nA,1), 'VariableNames', {'A', 'B'});
    validIdx = ~isnan(pairedIdxA);
    pairTable.B(validIdx) = B(pairedIdxA(validIdx));

    % Output struct
    result.pairedIdxA    = pairedIdxA;
    result.unmatchedA    = unmatchedA;
    result.unmatchedB     = unmatchedB;
    result.multiMatchedB = multiMatchedB;
    result.pairTable     = pairTable;
end