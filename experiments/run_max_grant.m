function out = run_max_grant(params)
% RUN_MAX_GRANT  Largest revenue-balanceable Grunderbe at the surtax cap.
%
%   out = run_max_grant([params])
%
% The no-exemption surtax tau_add is capped at 1 - tau0 (top marginal estate
% rate 100%). The largest fundable grant solves the fixed point
%
%     G* = tau_add_cap * E[b * death_rate] / N_entry     (in EUR)
%
% where the bequest flow itself depends on G through the stationary
% distribution. Iterates KFE-only solves at the cap (HJB is invariant to
% tau_add and G under the Step 1 gross-estate warm glow). Saves the grant
% steady state AND the status quo (m_sq) for distribution comparisons.

if nargin < 1 || isempty(params), params = params_default(); end
params.kfe_method = 'iterative';
cap = 1 - params.tau0;

grids = grids_build(params);
ops   = operator_build(grids, params);
[V, c, adrift, A, n_iter] = hjb_solve(grids, ops, params);
fprintf('HJB solved once (%d iters); fixed point on G at tau_add = %.2f.\n', n_iter, cap);

    function [m, rev, entry] = kfe_at(tadd, Geur)
        pp = params; pp.G_eur = Geur; pp.tau_add = tadd; pp = params_derive(pp);
        m = kfe_solve(grids, ops, pp, A);
        rev   = tadd * sum(grids.aa .* ops.death_rate .* m);
        entry = sum(ops.death_rate .* m);
    end

% Status quo (G = 0, no surtax) for comparison
[m_sq, ~, ~] = kfe_at(0, 0);

G = 2e5;   % start near the known frontier
for it = 1:20
    tk = tic;
    [m, rev, entry] = kfe_at(cap, G);
    G_new = rev / entry * params.eur_per_unit;
    fprintf('  [it %d] G = %8.0f -> balanceable G = %8.0f (%.0f s)\n', it, G, G_new, toc(tk));
    if abs(G_new - G) < 1e-3 * G_new, G = G_new; break; end
    G = G_new;
end

out = struct('G_max_eur', G, 'tau_add', cap, 'm', m, 'm_sq', m_sq, ...
             'params', params, 'hjb_iter', n_iter);
sol = struct('grids', grids, 'ops', ops, 'm', m, 'c', c, 'V', V, 'params', params);
out.mom    = moments(setfield(sol, 'm', m),    params); %#ok<SFLD>
out.mom_sq = moments(setfield(sol, 'm', m_sq), params); %#ok<SFLD>

fprintf('\n--- Largest balanceable Grunderbe (mortality=%s): G_max = %.0f EUR at tau_add = %.2f ---\n', ...
    params.mortality, G, cap);
fprintf('max grant : bottom 50%% %.3f, Gini %.3f, mean wealth %.0f EUR\n', ...
    out.mom.bottom50, out.mom.gini, out.mom.wealth_mean_eur);
fprintf('status quo: bottom 50%% %.3f, Gini %.3f, mean wealth %.0f EUR\n', ...
    out.mom_sq.bottom50, out.mom_sq.gini, out.mom_sq.wealth_mean_eur);

out_dir = fullfile(fileparts(fileparts(mfilename('fullpath'))), 'results', params.mortality);
if ~exist(out_dir, 'dir'), mkdir(out_dir); end
m_grant = m; %#ok<NASGU>
save(fullfile(out_dir, sprintf('maxgrant_%s.mat', params.mortality)), ...
     'params', 'm_grant', 'm_sq', 'out', '-v7');
fprintf('saved results to %s\n', out_dir);

end
