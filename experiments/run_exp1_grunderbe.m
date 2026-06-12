function [sol, mom, tau_star] = run_exp1_grunderbe(params, tau_lo, tau_hi)
% RUN_EXP1_GRUNDERBE  Experiment 1: Grunderbe with INCREMENTAL revenue balance.
%
%   [sol, mom, tau_star] = run_exp1_grunderbe()        % default params
%   [sol, mom, tau_star] = run_exp1_grunderbe(params)
%   [sol, mom, tau_star] = run_exp1_grunderbe(params, tau_lo, tau_hi)  % resume bracket
%
% G = 20,000 EUR at entry (h = 18), funded by raising the Step 1 flat rate
% tau0 until the revenue INCREMENT over the status quo pays the grant bill
% in the new stationary state:
%
%     Rev(tau*) = Rev_status_quo + G * N_entry.
%
% The status-quo revenue flow keeps leaking out of the household sector in
% BOTH scenarios (it implicitly funds unmodelled government spending), so
% the experiment isolates redistribution rather than mixing it with the
% recycling of previously-leaked revenue. NOTE: this deviates from the
% appendix §3.4 condition (total revenue = grant bill) by user decision
% 2026-06-12 — see CLAUDE.md discrepancy 13.
%
% Root-finding: Illinois method on tau (revenue is smooth and monotone).
% Tolerance: |increment - bill| < 0.1% of the bill (§5.3). HJB solves are
% warm-started. Saves results and comparison figures to results/step1/.

if nargin < 1 || isempty(params), params = params_default(); end
if nargin < 2, tau_lo = params.tau0; end   % status-quo rate: increment >= 0
if nargin < 3, tau_hi = 0.6; end

out_dir = fullfile(fileparts(fileparts(mfilename('fullpath'))), 'results', 'step1');
if ~exist(out_dir, 'dir'), mkdir(out_dir); end

% Status-quo revenue flow (the leak to be held constant across scenarios)
sq_file = fullfile(out_dir, 'status_quo.mat');
if exist(sq_file, 'file')
    sq = load(sq_file, 'mom');
    rev_sq = sq.mom.revenue;
    fprintf('status-quo revenue %.5f loaded from %s\n', rev_sq, sq_file);
else
    fprintf('status_quo.mat not found; solving the status quo first...\n');
    sol_sq = equilibrium(params);
    rev_sq = moments(sol_sq, params).revenue;
end

params.G_eur = 2e4;
params = params_derive(params);

t0 = tic;
V_warm = []; sol = []; mom = []; bill = NaN;

    function f = budget_resid(tau)
        params.tau0 = tau;
        sol = equilibrium(params, V_warm);
        V_warm = sol.V;
        mom  = moments(sol, params);
        bill = params.G * mom.entry_flow;
        f = mom.revenue - rev_sq - bill;
        fprintf('  tau0 = %.5f -> increment %.5f vs grant bill %.5f (budget %+.2e, hjb %d)\n', ...
            tau, mom.revenue - rev_sq, bill, f, sol.hjb_iter);
    end

f_lo = budget_resid(tau_lo);
tol  = 1e-3 * bill;
if f_lo > -tol  % grant already funded at the lower bracket edge
    tau_star = tau_lo;
else
    f_hi = budget_resid(tau_hi);
    if f_hi < 0
        error('run_exp1_grunderbe:bracket', 'tau_hi = %.3f cannot fund the grant', tau_hi);
    end
    side = 0; tau_star = NaN;
    for it = 1:30
        tau_new = tau_hi - f_hi * (tau_hi - tau_lo) / (f_hi - f_lo);
        f_new   = budget_resid(tau_new);
        tol     = 1e-3 * bill;
        if abs(f_new) < tol
            tau_star = tau_new;
            break
        end
        if f_new < 0
            tau_lo = tau_new; f_lo = f_new;
            if side == -1, f_hi = f_hi / 2; end
            side = -1;
        else
            tau_hi = tau_new; f_hi = f_new;
            if side == +1, f_lo = f_lo / 2; end
            side = +1;
        end
    end
    if isnan(tau_star)
        error('run_exp1_grunderbe:nobalance', 'incremental balance not reached in 30 iterations');
    end
end

fprintf('\n--- Experiment 1: Grunderbe, incremental balance (Step 1) ---\n');
fprintf('runtime %.0f s; tau* = %.4f (status-quo tau0 = %.2f, exemption %.0f EUR, grant %.0f EUR)\n', ...
    toc(t0), tau_star, params_default().tau0, params.F_eur, params.G_eur);
fprintf('top 1%% / 10%% %0.3f / %.3f, bottom 50%% %.3f, Gini %.3f, mean wealth %.0f EUR\n', ...
    mom.top1, mom.top10, mom.bottom50, mom.gini, mom.wealth_mean_eur);
if exist('sq', 'var')
    fprintf('status quo:   top 1%% / 10%% %.3f / %.3f, bottom 50%% %.3f, Gini %.3f, mean wealth %.0f EUR\n', ...
        sq.mom.top1, sq.mom.top10, sq.mom.bottom50, sq.mom.gini, sq.mom.wealth_mean_eur);
end

V = sol.V; m = sol.m; c = sol.c; %#ok<NASGU>
save(fullfile(out_dir, 'exp1_grunderbe.mat'), 'params', 'V', 'm', 'c', 'mom', 'tau_star', 'rev_sq');

% Comparison figures vs. status quo
if exist(sq_file, 'file')
    sq_full = load(sq_file, 'mom');
    g  = sol.grids;
    a_eur = g.a * params.eur_per_unit;

    f = figure('visible', 'off'); hold on;
    semilogy(a_eur(2:end)/1e3, sq_full.mom.wealth_marginal(2:end) ./ g.wa(2:end), 'LineWidth', 1.2);
    semilogy(a_eur(2:end)/1e3, mom.wealth_marginal(2:end) ./ g.wa(2:end), 'LineWidth', 1.2);
    set(gca, 'YScale', 'log'); legend('status quo', 'Grunderbe');
    xlabel('Wealth (1000 EUR)'); ylabel('density g(a)');
    title('Exp 1 (incremental balance): stationary wealth density');
    exportgraphics(f, fullfile(out_dir, 'exp1_wealth_density.png'), 'Resolution', 150);

    f = figure('visible', 'off'); hold on;
    plot(g.h, sq_full.mom.wealth_by_age * params.eur_per_unit / 1e3, 'LineWidth', 1.2);
    plot(g.h, mom.wealth_by_age * params.eur_per_unit / 1e3, 'LineWidth', 1.2);
    xline(params.hR, '--', 'HandleVisibility', 'off');
    legend('status quo', 'Grunderbe', 'Location', 'northwest');
    xlabel('Age'); ylabel('Mean wealth (1000 EUR)');
    title('Exp 1 (incremental balance): wealth-by-age profile');
    exportgraphics(f, fullfile(out_dir, 'exp1_wealth_by_age.png'), 'Resolution', 150);
end

fprintf('saved results and figures to %s\n', out_dir);

end
