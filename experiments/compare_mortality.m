function res = compare_mortality(params)
% COMPARE_MORTALITY  Step 1 (constant) vs Step 2.C (Destatis) mortality.
%
%   res = compare_mortality([params])
%
% Solves the status quo under both mortality specifications and reports the
% death-age distribution and bequest moments, to show how age-rising mortality
% fixes the inheritance flow (deaths concentrated at old ages instead of the
% spurious aging-out pile-up at h_max). Saves a comparison figure.

if nargin < 1, params = params_default(); end
out_dir = fullfile(fileparts(fileparts(mfilename('fullpath'))), 'results', 'step1');
if ~exist(out_dir, 'dir'), mkdir(out_dir); end

modes = {'step1', 'step2'};
labels = {'Step 1 (constant 2%)', 'Step 2.C (Destatis)'};
res = struct();

fprintf('\n%-26s %10s %10s\n', 'metric', labels{1}, labels{2});
for k = 1:2
    p = params; p.mortality = modes{k};
    sol = equilibrium(p);
    g = sol.grids; m = sol.m; ops = sol.ops; eur = p.eur_per_unit;
    dr = ops.death_rate; Nd = sum(dr .* m);

    R.deaths_by_age = accumarray(g.ih, dr .* m, [g.Nh, 1]);
    R.wealth_by_age = accumarray(g.ih, g.aa .* m, [g.Nh, 1]) ./ accumarray(g.ih, m, [g.Nh, 1]);
    R.mean_age_death = sum(g.hh .* dr .* m) / Nd;
    R.agingout_share = sum(ops.exit .* m) / Nd;
    R.mean_beq_eur   = sum(g.aa .* dr .* m) / Nd * eur;
    [asrt, ord] = sort(g.aa); cdf = cumsum(dr(ord) .* m(ord)) / Nd;
    R.med_beq_eur    = asrt(find(cdf >= 0.5, 1)) * eur;
    R.mean_wealth_eur = sum(g.aa .* m) * eur;
    R.bwr = sum(g.aa .* dr .* m) / sum(g.aa .* m);
    R.h = g.h; R.hjb_iter = sol.hjb_iter; R.kfe_res = sol.kfe.residual;
    res.(modes{k}) = R;
end

s1 = res.step1; s2 = res.step2;
row = @(name, a, b, f) fprintf('%-26s %10s %10s\n', name, sprintf(f,a), sprintf(f,b));
row('mean age at death',     s1.mean_age_death, s2.mean_age_death, '%.1f');
row('aging-out @h_max share', 100*s1.agingout_share, 100*s2.agingout_share, '%.1f%%');
row('mean bequest (EUR)',    s1.mean_beq_eur,   s2.mean_beq_eur,   '%.0f');
row('median bequest (EUR)',  s1.med_beq_eur,    s2.med_beq_eur,    '%.0f');
row('mean wealth (EUR)',     s1.mean_wealth_eur, s2.mean_wealth_eur, '%.0f');
row('bequest/wealth ratio',  s1.bwr,            s2.bwr,            '%.4f');
row('HJB iters',             s1.hjb_iter,       s2.hjb_iter,       '%.0f');

% Figure: death-age distribution and wealth-by-age, both specs
f = figure('visible', 'off');
subplot(1,2,1); hold on;
plot(s1.h, s1.deaths_by_age, 'LineWidth', 1.4);
plot(s2.h, s2.deaths_by_age, 'LineWidth', 1.4);
xline(params.h_max, ':', 'h_{max}');
xlabel('Age'); ylabel('death flow (mass/yr)'); legend(labels, 'Location', 'northwest');
title('Deaths by age');
subplot(1,2,2); hold on;
plot(s1.h, s1.wealth_by_age * params.eur_per_unit/1e3, 'LineWidth', 1.4);
plot(s2.h, s2.wealth_by_age * params.eur_per_unit/1e3, 'LineWidth', 1.4);
xline(params.hR, '--'); xlabel('Age'); ylabel('mean wealth (1000 EUR)');
legend(labels, 'Location', 'northwest'); title('Wealth by age');
exportgraphics(f, fullfile(out_dir, 'compare_mortality.png'), 'Resolution', 150);
close(f);
fprintf('saved figure: results/step1/compare_mortality.png\n');

end
