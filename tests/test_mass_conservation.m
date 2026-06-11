function test_mass_conservation()
% TEST_MASS_CONSERVATION  Boundary-operator checks on the full Step 1 model
% (appendix 5.1.2 and 5.1.6).
%
% Solves the full (a, y, h) stationary equilibrium on a reduced wealth grid
% (conservation properties are grid-size independent) and checks:
%   1. the generator-plus-kernel conserves mass exactly (column sums ~ 0),
%   2. the stationary solve has a small residual,
%   3. masses are nonnegative and sum to 1,
%   4. entry inflow at h0 equals the total death + aging-out outflow,
%   5. sanity: the wealth-by-age profile is hump-shaped (peak in the
%      interior, not at entry).

p = params_default();
p.Na = 100;   % reduced for runtime; checks are exact at any size

grids = grids_build(p);
ops   = operator_build(grids, p);

[~, ~, ~, A, n_iter] = hjb_solve(grids, ops, p);
fprintf('HJB converged in %d iterations\n', n_iter);

[m, ~, info] = kfe_solve(grids, ops, p, A);

% 1. Exact conservation of the discrete operator
scale = max(abs(ops.death_rate)) + 1 / grids.dh;
fprintf('max |column sum| of (A''+R): %.2e\n', info.colsum_max);
assert(info.colsum_max < 1e-10 * scale, 'generator does not conserve mass');

% 2. Stationarity residual
fprintf('stationary residual max|(A''+R) m|: %.2e\n', info.residual);
assert(info.residual < 1e-9, 'KFE residual too large');

% 3. Proper distribution
fprintf('sum(m) - 1 = %.2e, min(m) = %.2e\n', sum(m) - 1, min(m));
assert(abs(sum(m) - 1) < 1e-12, 'mass not normalised');
assert(min(m) > -1e-12, 'negative masses');

% 4. Entry inflow balances bequest-leaving outflow
inflow  = sum(info.R * m);
outflow = sum(ops.death_rate .* m);
fprintf('entry inflow %.6e vs death outflow %.6e (rel diff %.2e)\n', ...
    inflow, outflow, abs(inflow - outflow) / outflow);
assert(abs(inflow - outflow) < 1e-12 + 1e-10 * outflow, ...
    'inheritance kernel does not balance mortality outflow');

% 5. Wealth-by-age sanity (dev order item 3)
mass_h   = accumarray(grids.ih, m);
wealth_h = accumarray(grids.ih, grids.aa .* m) ./ mass_h;
[~, i_peak] = max(wealth_h);
fprintf('mean wealth by age: entry %.3f, peak %.3f at age %d, h_max %.3f\n', ...
    wealth_h(1), wealth_h(i_peak), grids.h(i_peak), wealth_h(end));
assert(i_peak > 1 && i_peak < grids.Nh, ...
    'wealth-by-age profile not hump-shaped');

disp('PASS test_mass_conservation');

end
