function out = run_exp1_flat_grant(params, G_eur)
% RUN_EXP1_FLAT_GRANT  Experiment 1 variant: grant funded by a no-exemption surtax.
%
%   out = run_exp1_flat_grant([params], [G_eur])   % default G = 100,000 EUR
%
% Funds the Grunderbe with a flat surtax on ALL bequests (no exemption), on top
% of the unchanged status-quo estate tax. Core solver: solve_grant_flat_tax.m.
% Reports feasibility and the required surtax rate; saves results and figures.

if nargin < 1 || isempty(params), params = params_default(); end
if nargin < 2 || isempty(G_eur),  G_eur  = 1e5; end

out_dir = fullfile(fileparts(fileparts(mfilename('fullpath'))), 'results', 'step1');
if ~exist(out_dir, 'dir'), mkdir(out_dir); end

glabel = sprintf('%gk EUR', G_eur/1e3);
tag    = matlab.lang.makeValidName(sprintf('flatgrant_%gk', G_eur/1e3));

t0  = tic;
out = solve_grant_flat_tax(params, G_eur);
mom = out.mom; sq = out.mom_sq;

fprintf('\n--- Grunderbe %s, NO-EXEMPTION surtax (mortality=%s) ---\n', glabel, params.mortality);
if out.feasible
    series_label = sprintf('Grunderbe %s (flat surtax)', glabel);
    fprintf('FEASIBLE: surtax tau_add* = %.4f on every bequest (%d KFE solves, %.0f s)\n', ...
        out.tau_add, out.n_kfe, toc(t0));
    fprintf('  surtax revenue %.5f = grant bill %.5f\n', out.surtax_rev, out.bill);
    fprintf('  grunderbe : top 1%% / 10%% %.3f / %.3f, bottom 50%% %.3f, Gini %.3f, mean wealth %.0f EUR\n', ...
        mom.top1, mom.top10, mom.bottom50, mom.gini, mom.wealth_mean_eur);
    fprintf('  status quo: top 1%% / 10%% %.3f / %.3f, bottom 50%% %.3f, Gini %.3f, mean wealth %.0f EUR\n', ...
        sq.top1, sq.top10, sq.bottom50, sq.gini, sq.wealth_mean_eur);
else
    series_label = sprintf('%s grant @ tau_add=%.2f (deficit)', glabel, out.tau_add_max);
    fprintf('INFEASIBLE even at tau_add=%.2f (%d KFE solves, %.0f s).\n', out.tau_add_max, out.n_kfe, toc(t0));
    fprintf('  surtax revenue %.5f < grant bill %.5f; largest fundable grant ~%.0f EUR.\n', ...
        out.surtax_rev, out.bill, out.max_grant_eur);
end

sol = out.sol; V = sol.V; m = sol.m; c = sol.c; %#ok<NASGU>
save(fullfile(out_dir, [tag '.mat']), 'params', 'G_eur', 'V', 'm', 'c', 'mom', 'out');
save_experiment_figures(sol, mom, sq, out_dir, tag, series_label, ...
    sprintf('Exp 1 flat surtax (Grunderbe %s)', glabel));
fprintf('saved results and figures to %s\n', out_dir);

end
