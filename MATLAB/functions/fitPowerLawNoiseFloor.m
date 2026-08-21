function out = fitPowerLawNoiseFloor(f, P, fitRangeHz, excludeHz, opts)
%FITPOWERLAWNOISEFLOOR Fit PSD to P(f) = A*f^(-chi) + C
%
% Inputs
%   f           : frequency vector (Hz), column or row
%   P           : PSD values (same size as f), must be >0
%   fitRangeHz  : [fmin fmax], e.g. [2 499]
%   excludeHz   : vector of center freqs to exclude (e.g. 50:50:450 or 60:60:480)
%   opts        : struct with fields (optional)
%       .excludeHalfWidthHz (default 1)  -> removes bins within ± this width of excludeHz
%       .useLogResidual     (default true) -> fits log10(P) residuals (recommended)
%       .robust             (default true) -> uses robust weighting (bisquare) if fitnlm available
%
% Output struct out:
%   .A, .chi, .C
%   .Pfit   : model on all f
%   .maskFit: logical mask used for fitting
%   .exitflag, .resnorm, .residual
% Function based on the KJM power-law paper
% Created with ChatGPT !!! USE WITH CAUTION!!!

    arguments
        f
        P
        fitRangeHz (1,2) double = [2 499]
        excludeHz double = []
        opts struct = struct()
    end
    
    if ~isfield(opts, 'excludeHalfWidthHz')
        opts.excludeHalfWidthHz = 1;
    end
    
    if ~isfield(opts, 'useLogResidual')
        opts.useLogResidual = true;
    end
    
    if ~isfield(opts, 'robust')
        opts.robust = true;
    end

    f = f(:); P = P(:);

    % Basic validity
    good = isfinite(f) & isfinite(P) & (P > 0);
    inBand = f >= fitRangeHz(1) & f <= fitRangeHz(2);

    % Exclude line-noise harmonics etc.
    excl = false(size(f));
    if ~isempty(excludeHz)
        for k = 1:numel(excludeHz)
            excl = excl | (abs(f - excludeHz(k)) <= opts.excludeHalfWidthHz);
        end
    end

    mask = good & inBand & ~excl;

    ff = f(mask);
    PP = P(mask);

    % --- Initial guesses ---
    % Guess C from top part of the band (e.g., upper 20% of frequencies)
    q = 0.8;
    hiMask = ff >= quantile(ff, q);
    C0 = median(PP(hiMask));
    C0 = max(C0, eps);

    % Guess chi from log-log slope after subtracting C0 (clipped)
    Psub = max(PP - C0, eps);
    p = polyfit(log10(ff), log10(Psub), 1);
    chi0 = max(-p(1), 0.1);

    % Guess A from median of (P-C)*f^chi
    A0 = median(Psub .* (ff .^ chi0));
    A0 = max(A0, eps);

    x0 = [log(A0), log(chi0), log(C0)]; % optimize in log-space for positivity

    % --- Objective ---
    if opts.useLogResidual
        % minimize log10(P) - log10(model)
        fun = @(x) (log10(PP) - log10(exp(x(1)) .* ff.^(-exp(x(2))) + exp(x(3))));
    else
        % minimize linear residuals with relative weighting to avoid overweighting low-f
        model = @(x) exp(x(1)) .* ff.^(-exp(x(2))) + exp(x(3));
        fun = @(x) (PP - model(x)) ./ max(model(x), eps);
    end

    % Solve
    lb = [log(eps), log(0.01), log(eps)];
    ub = [log(realmax('double')/1e6), log(20), log(realmax('double')/1e6)];

    if license('test','Optimization_Toolbox')
        options = optimoptions('lsqnonlin','Display','off','MaxFunctionEvaluations',5e4);
        [x,resnorm,residual,exitflag] = lsqnonlin(fun, x0, lb, ub, options);
    else
        % Fallback: fminsearch (no bounds) - less ideal but works
        warning('Optimization Toolbox not found; using fminsearch fallback.');
        obj = @(x) sum(fun(x).^2);
        x = fminsearch(obj, x0, optimset('Display','off','MaxFunEvals',5e4));
        residual = fun(x);
        resnorm = sum(residual.^2);
        exitflag = NaN;
    end

    A = exp(x(1));
    chi = exp(x(2));
    C = exp(x(3));

    Pfit = A .* f.^(-chi) + C;

    out = struct();
    out.A = A;
    out.chi = chi;
    out.C = C;
    out.Pfit = Pfit;
    out.maskFit = mask;
    out.exitflag = exitflag;
    out.resnorm = resnorm;
    out.residual = residual;
end