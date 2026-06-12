function [sol, mom, tau_star] = run_exp1_grunderbe(params, tau_lo, tau_hi)
% RUN_EXP1_GRUNDERBE  Experiment 1: Grunderbe with revenue balance (§3.4, outline 6.1).
%
%   [sol, mom, tau_star] = run_exp1_grunderbe()        % default params
%   [sol, mom, tau_star] = run_exp1_grunderbe(params)
%   [sol, mom, tau_star] = run_exp1_grunderbe(params, tau_lo, tau_hi)  % resume bracket
%
% G = 20,000 EUR at entry (h = 18), funded by bisecting on the Step 1 flat
% estate-tax rate tau0 until revenue = G * N_entry in the NEW stationary
% state (appendix §3.4; tolerance 0.1% of the grant bill, §5.3). Inner loop
% is a full HJB+KFE solve. Saves results and comparison figures vs. the
% status quo (if results/step1/status_quo.mat exists) to results/step1/.

if nargin < 1 || isempty(params), params = params_default(); end
if nargin < 2, tau_lo = 0; end
if nargin < 3, tau_hi = 0.6; end
params.G_eur = 2e4;
params = params_derive(params);

out_dir = fullfile(fileparts(fileparts(mfilename('fullpath'))), 'results', 'step1');
if ~exist(out_dir, 'dir'), mkdir(out_dir); end

sol = []; mom = []; tau_star = NaN;
t0 = tic;
V_warm = [];
for outer = 1:50
    tau_mid = 0.5 * (tau_lo + tau_hi);
    params.tau0 = tau_mid;
    sol = equilibrium(params, V_warm);
    V_warm = sol.V;
    mom = moments(sol, params);
    bill   = params.G * mom.entry_flow;
    budget = mom.revenue - bill;
    fprintf('  tau0 = %.5f -> revenue %.5f vs grant bill %.5f (budget %+.2e)\n', ...
        tau_mid, mom.revenue, bill, budget);
    if abs(budget) < 1e-3 * bill
        tau_star = tau_mid;
        break
    end
    if budget < 0, tau_lo = tau_mid; else, tau_hi = tau_mid; end
end
if isnan(tau_star)
    error('run_exp1_grunderbe:nobalance', 'revenue balance not reached in 50 bisections');
end

fprintf('\n--- Experiment 1: Grunderbe (Step 1) ---\n');
fprintf('runtime %.0f s; revenue-balancing tau0 = %.4f (exemption %.0f EUR, grant %.0f EUR)\n', ...
    toc(t0), tau_star, params.F_eur, params.G_eur);
fprintf('top 1%% / 10%% %0.3f / %.3f, bottom 50%% %.3f, Gini %.3f, mean wealth %.0f EUR\n', ...
    mom.top1, mom.top10, mom.bottom50, mom.gini, mom.wealth_mean_eur);

V = sol.V; m = sol.m; c = sol.c; %#ok<NASGU>
save(fullfile(out_dir, 'exp1_grunderbe.mat'), 'params', 'V', 'm', 'c', 'mom', 'tau_star');

% Comparison figures vs. status quo
sq_file = fullfile(out_dir, 'status_quo.mat');
if exist(sq_file, 'file')
    sq = load(sq_file, 'mom', 'params');
    g  = sol.grids;
    a_eur = g.a * params.eur_per_unit;

    f = figure('visible', 'off'); hold on;
    semilogy(a_eur(2:end)/1e3, sq.mom.wealth_marginal(2:end) ./ g.wa(2:end), 'LineWidth', 1.2);
    semilogy(a_eur(2:end)/1e3, mom.wealth_marginal(2:end) ./ g.wa(2:end), 'LineWidth', 1.2);
    set(gca, 'YScale', 'log'); legend('status quo', 'Grunderbe');
    xlabel('Wealth (1000 EUR)'); ylabel('density g(a)');
    title('Exp 1: stationary wealth density');
    exportgraphics(f, fullfile(out_dir, 'exp1_wealth_density.png'), 'Resolution', 150);

    f = figure('visible', 'off'); hold on;
    plot(g.h, sq.mom.wealth_by_age * params.eur_per_unit / 1e3, 'LineWidth', 1.2);
    plot(g.h, mom.wealth_by_age * params.eur_per_unit / 1e3, 'LineWidth', 1.2);
    xline(params.hR, '--', 'HandleVisibility', 'off');
    legend('status quo', 'Grunderbe', 'Location', 'northwest');
    xlabel('Age'); ylabel('Mean wealth (1000 EUR)');
    title('Exp 1: wealth-by-age profile');
    exportgraphics(f, fullfile(out_dir, 'exp1_wealth_by_age.png'), 'Resolution', 150);
end

fprintf('saved results and figures to %s\n', out_dir);

end
