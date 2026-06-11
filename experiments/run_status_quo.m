function [sol, mom] = run_status_quo(params)
% RUN_STATUS_QUO  Step 1 status-quo stationary equilibrium (dev order item 5).
%
%   [sol, mom] = run_status_quo()          % full-resolution default params
%   [sol, mom] = run_status_quo(params)    % custom params
%
% Solves the stationary equilibrium under the current (placeholder) estate
% tax with no grant, prints summary moments, and saves results + figures to
% results/step1/.

if nargin < 1, params = params_default(); end

out_dir = fullfile(fileparts(fileparts(mfilename('fullpath'))), 'results', 'step1');
if ~exist(out_dir, 'dir'), mkdir(out_dir); end

t0 = tic;
sol = equilibrium(params);
mom = moments(sol, params);
runtime = toc(t0);

fprintf('\n--- Status quo (Step 1) ---\n');
fprintf('runtime %.1f s (HJB %d iterations), KFE residual %.1e\n', ...
    runtime, sol.hjb_iter, sol.kfe.residual);
fprintf('mean wealth        %8.2f units = %.0f EUR\n', ...
    mom.wealth_mean, mom.wealth_mean_eur);
fprintf('top 1%% / 5%% / 10%%  %6.3f / %.3f / %.3f\n', mom.top1, mom.top5, mom.top10);
fprintf('bottom 50%%         %8.3f\n', mom.bottom50);
fprintf('Gini               %8.3f\n', mom.gini);
fprintf('bequest/wealth     %8.4f per year\n', mom.bwr);
fprintf('tax revenue        %8.4f units per year (%.2f%% of bequest flow)\n', ...
    mom.revenue, 100 * mom.revenue / (mom.bwr * mom.wealth_mean));

% Persist (regenerable pieces dropped: A, ops, grids)
V = sol.V; m = sol.m; c = sol.c; %#ok<NASGU>
save(fullfile(out_dir, 'status_quo.mat'), 'params', 'V', 'm', 'c', 'mom');

% Figures
g = sol.grids;

f = figure('visible', 'off');
ga_density = mom.wealth_marginal ./ g.wa;       % marginal density on a
a_eur = g.a * params.eur_per_unit;
semilogy(a_eur(2:end) / 1e3, ga_density(2:end), 'LineWidth', 1.5);
xlabel('Wealth (1000 EUR)'); ylabel('density g(a)');
title('Status quo: stationary wealth density');
exportgraphics(f, fullfile(out_dir, 'status_quo_wealth_density.png'), 'Resolution', 150);

f = figure('visible', 'off');
plot(g.h, mom.wealth_by_age * params.eur_per_unit / 1e3, 'LineWidth', 1.5);
xline(params.hR, '--'); xlabel('Age'); ylabel('Mean wealth (1000 EUR)');
title('Status quo: wealth-by-age profile');
exportgraphics(f, fullfile(out_dir, 'status_quo_wealth_by_age.png'), 'Resolution', 150);

fprintf('saved results and figures to %s\n', out_dir);

end
