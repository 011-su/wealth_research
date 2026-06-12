function sol = equilibrium(params, V0)
% EQUILIBRIUM  Stationary equilibrium for given params (appendix 3.2-3.3).
%
%   sol = equilibrium(params)
%   sol = equilibrium(params, V0)   % warm-start the HJB (same grid sizes)
%
% Wrapper: grids -> operators -> HJB -> KFE. Returns a struct with grids,
% ops, V, c, adrift, A, masses m, density g, KFE diagnostics, and the
% estate-tax revenue flow (per unit time, model units).
%
% The revenue-balance bisection on tau (appendix 3.4) is NOT here yet; it
% arrives with run_exp1_grunderbe.m (development order item 7).

if nargin < 2, V0 = []; end
grids = grids_build(params);
ops   = operator_build(grids, params);

[V, c, adrift, A, n_iter] = hjb_solve(grids, ops, params, V0);
[m, g, kfe_info]          = kfe_solve(grids, ops, params, A);

sol = struct('grids', grids, 'ops', ops, 'V', V, 'c', c, ...
             'adrift', adrift, 'A', A, 'm', m, 'g', g, ...
             'hjb_iter', n_iter, 'kfe', kfe_info);

% Estate-tax revenue: E[T_e(a) * death rate] under the stationary masses
sol.revenue = sum(ops.tax_fn(grids.aa, params) .* ops.death_rate .* m);

end
