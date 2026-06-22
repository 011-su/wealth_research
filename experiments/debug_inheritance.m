function info = debug_inheritance(params)
% DEBUG_INHERITANCE  Decompose the bequest / inheritance flow in the status quo.
%
%   info = debug_inheritance()         % uses cached results/step1/status_quo.mat
%   info = debug_inheritance(params)   % solves fresh if no cache / custom params
%
% Built to answer "why does the inheritance flow / estate-tax base look small?"
% Reports, in model units AND EUR, every relevant cut of the bequest flow:
%   - bequest LEFT per decedent (gross), split by death type (mortality vs the
%     forced aging-out at h_max);
%   - inheritance RECEIVED per heir (post-tax = starting wealth);
%   - the taxable slice above the exemption, and the effective tax rate;
%   - an aggregate scaling to Germany vs published benchmarks;
%   - the wealth-distribution context (skew, top shares) that thins the base.
% Saves a figure of the bequest-flow CDF over wealth with the exemption marked.

if nargin < 1, params = params_default(); end

proj    = fileparts(fileparts(mfilename('fullpath')));
out_dir = fullfile(proj, 'results', 'step1');
cache   = fullfile(out_dir, 'status_quo.mat');

if exist(cache, 'file')
    S = load(cache, 'params', 'm', 'c');
    params = S.params; m = S.m; c = S.c;
    grids  = grids_build(params);
    ops    = operator_build(grids, params);
    fprintf('loaded cached status quo (Na=%d, theta_b=%.3g)\n', params.Na, params.theta_b);
else
    fprintf('no cache; solving status quo (Na=%d)...\n', params.Na);
    sol = equilibrium(params); grids = sol.grids; ops = sol.ops; m = sol.m; c = sol.c;
end
sol = struct('grids', grids, 'ops', ops, 'm', m, 'c', c);

aa  = grids.aa;  eur = params.eur_per_unit;
dr  = ops.death_rate;  lam = ops.lambda;  ex = ops.exit;   % rates by state
F   = params.F;

% ---- Flows (per capita = per unit population mass, per year) ----
Nd      = sum(dr  .* m);   Nd_mort = sum(lam .* m);   Nd_age = sum(ex .* m);
Bgross  = sum(aa .* dr  .* m);
Bg_mort = sum(aa .* lam .* m);
Bg_age  = sum(aa .* ex  .* m);
W       = sum(aa .* m);                       % mean wealth (mass 1)
taxable = sum(max(aa - F, 0) .* dr .* m);     % flow above the exemption
rev     = sum(ops.tax_fn(aa, params) .* dr .* m);
frac_above = sum((aa > F) .* dr .* m) / Nd;   % share of decedents above F

% ---- Per-decedent and per-heir ----
mean_beq      = Bgross  / Nd;
mean_beq_mort = Bg_mort / max(Nd_mort, realmin);
mean_beq_age  = Bg_age  / max(Nd_age , realmin);

% death-weighted median bequest, and population median wealth
[asrt, ord] = sort(aa);
cdf_death = cumsum(dr(ord) .* m(ord)) / Nd;          med_beq    = asrt(find(cdf_death >= 0.5, 1));
cdf_pop   = cumsum(m(ord));                           med_wealth = asrt(find(cdf_pop   >= 0.5, 1));

% inheritance received per heir = mean starting wealth at h = h0 (post-tax, +G)
ih0 = grids.ih == 1;  mass_h0 = sum(m(ih0));
mean_entry = sum(aa(ih0) .* m(ih0)) / mass_h0;

% ---- Aggregate scaling to Germany (illustrative) ----
n_hh = 40.9e6;   % German households (Destatis); model mass 1 = population
agg_beq_bn = Bgross * eur * n_hh / 1e9;
agg_rev_bn = rev    * eur * n_hh / 1e9;
deaths_M   = Nd * n_hh / 1e6;

% ---- Report ----
e = @(x) x * eur;
fprintf('\n================ INHERITANCE-FLOW DIAGNOSTIC (status quo) ================\n');
fprintf('units: 1 model unit = %.0f EUR;  exemption F = %.0f EUR = %.2f units\n', eur, params.F_eur, F);

fprintf('\n-- Bequest LEFT per decedent --\n');
fprintf('  decedents/yr (per capita)      : %.4f  (avg generation %.0f yr)\n', Nd, 1/Nd);
fprintf('  MEAN bequest / decedent (gross): %8.2f u = %9.0f EUR\n', mean_beq, e(mean_beq));
fprintf('  MEDIAN bequest / decedent      : %8.2f u = %9.0f EUR\n', med_beq, e(med_beq));
fprintf('  mean wealth (whole population)  : %8.2f u = %9.0f EUR\n', W, e(W));
fprintf('  -> data benchmark: ~200k-400k EUR per inheritance\n');

fprintf('\n-- Death type: constant hazard lets many reach h_max=%d (low wealth) --\n', params.h_max);
fprintf('  mortality deaths : %.4f/yr (%4.1f%%), mean bequest %9.0f EUR\n', ...
    Nd_mort, 100*Nd_mort/Nd, e(mean_beq_mort));
fprintf('  aging-out at h_max: %.4f/yr (%4.1f%%), mean bequest %9.0f EUR  <- drags the average down\n', ...
    Nd_age, 100*Nd_age/Nd, e(mean_beq_age));

fprintf('\n-- Inheritance RECEIVED per heir (post-tax = starting wealth) --\n');
fprintf('  mean starting wealth at h0      : %8.2f u = %9.0f EUR\n', mean_entry, e(mean_entry));
fprintf('  (entry mass %.4f vs decedents %.4f: %.1e -- mass conservation)\n', ...
    mass_h0, Nd, abs(mass_h0 - Nd));

fprintf('\n-- The exemption thins the tax base --\n');
fprintf('  decedents with estate > F (400k): %.1f%%  (F sits at the %.0fth wealth pctile)\n', ...
    100*frac_above, 100*(1 - frac_above));
fprintf('  gross bequest flow  : %.4f u/yr = %8.0f EUR/capita/yr\n', Bgross, e(Bgross));
fprintf('  taxable flow > F     : %.4f u/yr = %8.0f EUR/capita/yr (%.1f%% of gross)\n', ...
    taxable, e(taxable), 100*taxable/Bgross);
fprintf('  estate-tax revenue   : %.4f u/yr = %8.0f EUR/capita/yr (%.1f%% of gross)\n', ...
    rev, e(rev), 100*rev/Bgross);

fprintf('\n-- Aggregate scaling (x %.1fM households; illustrative) --\n', n_hh/1e6);
fprintf('  aggregate bequest flow : %6.1f bn EUR/yr   (data ~200-400 bn/yr, Tiefensee-Grabka)\n', agg_beq_bn);
fprintf('  aggregate estate-tax   : %6.1f bn EUR/yr   (Destatis ErbSt ~8-11 bn/yr)\n', agg_rev_bn);
fprintf('  implied deaths/yr      : %6.2f M           (Germany ~1.05 M/yr)\n', deaths_M);

mom = moments(sol, params);
fprintf('\n-- Wealth distribution (why the slice above 400k is thin) --\n');
fprintf('  median/mean wealth ratio : %.2f  (lower => more right-skew)\n', med_wealth / W);
fprintf('  top 1%% / 10%% / Gini      : %.3f / %.3f / %.3f\n', mom.top1, mom.top10, mom.gini);
fprintf('  -> data: top1 ~0.27, top10 ~0.60, Gini ~0.77. Model tail too THIN =>\n');
fprintf('     the taxable base above 400k is understated (known Step 1 limitation).\n');
fprintf('==========================================================================\n');

% ---- Figure: cumulative bequest-flow share vs wealth, exemption marked ----
flow_a = accumarray(grids.ia, aa .* dr .* m, [grids.Na, 1]);
cum_share = cumsum(flow_a) / sum(flow_a);
share_below_F = interp1(grids.a, cum_share, F, 'linear', 1);
a_keur = grids.a * eur / 1e3;

f = figure('visible', 'off'); hold on;
plot(a_keur, cum_share, 'LineWidth', 1.6);
xline(params.F_eur/1e3, '--', sprintf('exemption %.0fk', params.F_eur/1e3));
yline(share_below_F, ':');
xlim([0, 2000]);
xlabel('Estate size (1000 EUR)');
ylabel('cumulative share of bequest flow');
title(sprintf('Bequest flow by estate size — %.0f%% is below the 400k exemption', 100*share_below_F));
exportgraphics(f, fullfile(out_dir, 'debug_inheritance_flow_cdf.png'), 'Resolution', 150);
close(f);
fprintf('saved figure: results/step1/debug_inheritance_flow_cdf.png\n');

info = struct('Nd', Nd, 'Bgross', Bgross, 'taxable', taxable, 'rev', rev, 'W', W, ...
    'mean_beq', mean_beq, 'mean_beq_eur', e(mean_beq), 'med_beq_eur', e(med_beq), ...
    'mean_entry_eur', e(mean_entry), 'frac_above', frac_above, ...
    'Nd_age_share', Nd_age/Nd, 'agg_beq_bn', agg_beq_bn, 'agg_rev_bn', agg_rev_bn, ...
    'share_below_F', share_below_F, 'top1', mom.top1, 'gini', mom.gini);

end
