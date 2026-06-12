function [theta_b, sol, mom, hist] = calibrate_step1(params, target_bwr, bracket)
% CALIBRATE_STEP1  Calibrate theta_b to the bequest-to-wealth ratio (§4.1).
%
%   [theta_b, sol, mom, hist] = calibrate_step1(params)
%   [theta_b, sol, mom, hist] = calibrate_step1(params, target_bwr)
%   [theta_b, sol, mom, hist] = calibrate_step1(params, target_bwr, bracket)
%
% Root-finds bwr(theta_b) = target_bwr with the Illinois (modified regula
% falsi) method on log10(theta_b) -- much faster than bisection here
% because bwr is smooth, monotone, and nearly flat in theta_b under the
% constant Step 1 hazard. Successive HJB solves are warm-started with the
% previous V. Inner loop = full HJB+KFE at the resolution in params; pass a
% reduced-Na params for speed and re-solve at full resolution afterwards.
%
% bracket (optional): struct with fields lg_lo, lg_hi, and optionally
% bwr_lo, bwr_hi from a previous run, to skip the endpoint evaluations.
%
% Stops when |bwr - target| < 1% of target. hist records (theta_b, bwr).

if nargin < 2 || isempty(target_bwr), target_bwr = params.target_bwr; end
if nargin < 3, bracket = struct(); end
if ~isfield(bracket, 'lg_lo'), bracket.lg_lo = -3; end
if ~isfield(bracket, 'lg_hi'), bracket.lg_hi =  3; end

hist = [];
V_warm = [];

lg_lo = bracket.lg_lo;  lg_hi = bracket.lg_hi;
if isfield(bracket, 'bwr_lo'), f_lo = bracket.bwr_lo - target_bwr;
else, f_lo = eval_resid(lg_lo); end
if isfield(bracket, 'bwr_hi'), f_hi = bracket.bwr_hi - target_bwr;
else, f_hi = eval_resid(lg_hi); end

if ~(f_lo < 0 && f_hi > 0)
    error('calibrate_step1:bracket', ...
        'target bwr %.4f outside bracket [%.4f, %.4f]', ...
        target_bwr, f_lo + target_bwr, f_hi + target_bwr);
end

% Illinois method: secant step, halve the retained endpoint's residual
% whenever the same side is kept twice (guarantees convergence)
side = 0;
for it = 1:30
    lg_new = lg_hi - f_hi * (lg_hi - lg_lo) / (f_hi - f_lo);
    f_new  = eval_resid(lg_new);
    if abs(f_new) < 0.01 * target_bwr
        break
    end
    if f_new < 0
        lg_lo = lg_new; f_lo = f_new;
        if side == -1, f_hi = f_hi / 2; end
        side = -1;
    else
        lg_hi = lg_new; f_hi = f_new;
        if side == +1, f_lo = f_lo / 2; end
        side = +1;
    end
end

theta_b = 10^lg_new;
params.theta_b = theta_b;
sol = equilibrium(params, V_warm);
mom = moments(sol, params);

    function f = eval_resid(lg_theta)
        params.theta_b = 10^lg_theta;
        s = equilibrium(params, V_warm);
        V_warm = s.V;
        mo = moments(s, params);
        f = mo.bwr - target_bwr;
        hist(end+1, :) = [params.theta_b, mo.bwr]; %#ok<AGROW>
        fprintf('  theta_b = %9.4g  ->  bwr = %.4f (hjb %d iters)\n', ...
            params.theta_b, mo.bwr, s.hjb_iter);
    end

end
