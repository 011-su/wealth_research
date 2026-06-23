function out = run_grunderbe(params, G_eur)
% RUN_GRUNDERBE  Solve, report, and plot a Grunderbe experiment for any grant.
%
%   out = run_grunderbe(params, G_eur)
%
% Generic engine behind the named Grunderbe experiments (run_exp1_grunderbe =
% 20k, run_exp1b_grunderbe_200k = 200k). The grant G_eur is the ONLY policy
% input; every label, message, filename, and reported number is derived from
% it and from the solver output, so changing the grant (or any param) needs
% no other edits. Incremental revenue balance and the one-HJB speedup live in
% solve_revenue_balance.m.
%
% Returns the full struct from solve_revenue_balance (fields .feasible,
% .tau_star, .sol, .mom, .sol_sq, .mom_sq, .increment, .bill, ...).

if nargin < 1 || isempty(params), params = params_default(); end
params.kfe_method = 'iterative';   % gmres + ilu; matches direct to machine
                                   % precision (validated Na=300), faster + low
                                   % memory on the costly KFE sweeps. Override
                                   % with params.kfe_method='direct' if needed.

out_dir = fullfile(fileparts(fileparts(mfilename('fullpath'))), 'results', 'step1');
if ~exist(out_dir, 'dir'), mkdir(out_dir); end

% Labels and filenames derived from the grant (no hardcoded magnitudes)
grant_label = format_eur(G_eur);
tag         = matlab.lang.makeValidName(sprintf('grunderbe_%s', grant_label));
ttl         = sprintf('Grunderbe %s (incremental balance)', grant_label);

t0  = tic;
out = solve_revenue_balance(params, G_eur, 'tau_max', 1.0);
mom = out.mom; sq = out.mom_sq;

fprintf('\n--- Grunderbe %s, incremental balance (Step 1) ---\n', grant_label);
if out.feasible
    series_label = sprintf('Grunderbe %s', grant_label);
    fprintf('runtime %.0f s (%d KFE solves); tau* = %.4f (status-quo tau0 = %.2f, exemption %s)\n', ...
        toc(t0), out.n_kfe, out.tau_star, params.tau0, format_eur(params.F_eur));
    fprintf('grunderbe:  top 1%% / 10%% %.3f / %.3f, bottom 50%% %.3f, Gini %.3f, mean wealth %.0f EUR\n', ...
        mom.top1, mom.top10, mom.bottom50, mom.gini, mom.wealth_mean_eur);
    fprintf('status quo: top 1%% / 10%% %.3f / %.3f, bottom 50%% %.3f, Gini %.3f, mean wealth %.0f EUR\n', ...
        sq.top1, sq.top10, sq.bottom50, sq.gini, sq.wealth_mean_eur);
else
    series_label = sprintf('%s grant @ tau=%.2f (deficit)', grant_label, out.tau_max);
    fprintf('INFEASIBLE under revenue balance (%d KFE solves, %.0f s).\n', out.n_kfe, toc(t0));
    fprintf('  grant bill %.5f vs max increment %.5f at tau = %.2f.\n', ...
        out.bill, out.increment, out.tau_max);
    fprintf('  A %s grant to every entrant cannot be funded by an estate tax with a\n', grant_label);
    fprintf('  %s exemption: the largest revenue-balanceable grant is ~%s.\n', ...
        format_eur(params.F_eur), format_eur(out.max_grant_eur));
    fprintf('  Figures show the tau = %.2f equilibrium with the full grant (NOT\n', out.tau_max);
    fprintf('  revenue-balanced; the unfunded part is implicitly deficit-financed).\n');
end

sol = out.sol; V = sol.V; m = sol.m; c = sol.c;
save(fullfile(out_dir, [tag '.mat']), 'params', 'G_eur', 'V', 'm', 'c', 'mom', 'out');
save_experiment_figures(sol, mom, sq, out_dir, tag, series_label, ttl);
fprintf('saved results and figures to %s\n', out_dir);

end


function s = format_eur(x)
if x >= 1e6
    s = sprintf('%gM EUR', x / 1e6);
elseif x >= 1e3
    s = sprintf('%gk EUR', x / 1e3);
else
    s = sprintf('%g EUR', x);
end
end
