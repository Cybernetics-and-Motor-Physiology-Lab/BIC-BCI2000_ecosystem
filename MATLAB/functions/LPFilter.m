function S_out = LPFilter(S_in, T)
% LPFilter implements the BCI2000-style exponential low-pass filter
% Inputs:
%   S_in : input signal (vector or matrix; time × channels)
%   T    : time constant (scalar > 0)
% Output:
%   S_out : filtered output, same size as S_in

    if T <= 0
        error('Time constant T must be positive.');
    end

    alpha = exp(-1 / T);
    beta = 1 - alpha;

    % Initialize output
    S_out = zeros(size(S_in));
    
    % First sample
    S_out(1,:) = beta * S_in(1,:);

    % Recursive filter
    for t = 2:size(S_in,1)
        S_out(t,:) = alpha * S_out(t-1,:) + beta * S_in(t,:);
    end
end
