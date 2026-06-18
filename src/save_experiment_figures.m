function save_experiment_figures(sol, mom, sq_mom, out_dir, stem, series_label, ttl)
% SAVE_EXPERIMENT_FIGURES  Wealth-density and wealth-by-age comparison plots.
%
%   save_experiment_figures(sol, mom, sq_mom, out_dir, stem, series_label, ttl)
%
% Writes <stem>_wealth_density.png and <stem>_wealth_by_age.png to out_dir,
% overlaying the experiment moments (mom) on the status quo (sq_mom).

g = sol.grids; p = sol.params;
a_eur = g.a * p.eur_per_unit;

f = figure('visible', 'off'); hold on;
semilogy(a_eur(2:end)/1e3, sq_mom.wealth_marginal(2:end) ./ g.wa(2:end), 'LineWidth', 1.2);
semilogy(a_eur(2:end)/1e3, mom.wealth_marginal(2:end) ./ g.wa(2:end), 'LineWidth', 1.2);
set(gca, 'YScale', 'log'); legend('status quo', series_label);
xlabel('Wealth (1000 EUR)'); ylabel('density g(a)');
title([ttl ': stationary wealth density']);
exportgraphics(f, fullfile(out_dir, [stem '_wealth_density.png']), 'Resolution', 150);
close(f);

f = figure('visible', 'off'); hold on;
plot(g.h, sq_mom.wealth_by_age * p.eur_per_unit / 1e3, 'LineWidth', 1.2);
plot(g.h, mom.wealth_by_age * p.eur_per_unit / 1e3, 'LineWidth', 1.2);
xline(p.hR, '--', 'HandleVisibility', 'off');
legend('status quo', series_label, 'Location', 'northwest');
xlabel('Age'); ylabel('Mean wealth (1000 EUR)');
title([ttl ': wealth-by-age profile']);
exportgraphics(f, fullfile(out_dir, [stem '_wealth_by_age.png']), 'Resolution', 150);
close(f);

end
