function mom = moments(sol, params)
% MOMENTS  Distributional statistics of a stationary solution (appendix §6).
%
%   mom = moments(sol, params)
%
% Input:  sol from equilibrium(): uses grids, ops, masses m.
% Output: mom struct
%   .wealth_mean        mean wealth (model units; *_eur fields in EUR)
%   .top1, .top5, .top10, .bottom50   wealth shares
%   .gini               wealth Gini
%   .bwr                bequest-to-wealth ratio per year:
%                       E[a * death_rate] / E[a]  (deaths only, no gifts)
%   .revenue            estate-tax revenue flow (model units per year)
%   .entry_flow         population entry flow (per year)
%   .wealth_by_age      (Nh x 1) mean wealth profile
%   .wealth_marginal    (Na x 1) masses on the wealth grid
%
% Top shares and Gini treat the wealth marginal as an atomic distribution
% on the a grid; the threshold point mass is split pro rata for shares.

g  = sol.grids;
m  = sol.m;

% Wealth marginal (masses by grid point, ascending in a)
ma = accumarray(g.ia, m, [g.Na, 1]);
wealth = ma .* g.a;
W = sum(wealth);

mom.wealth_mean     = W;                % total mass is 1
mom.wealth_mean_eur = W * params.eur_per_unit;

% Lorenz-based statistics
cum_pop    = cumsum(ma);
cum_wealth = cumsum(wealth);

    function s = top_share(p)
        % wealth share of the richest fraction p, splitting the threshold atom
        thresh = 1 - p;
        k = find(cum_pop >= thresh, 1);
        below = 0; if k > 1, below = cum_wealth(k - 1); end
        inside = (cum_pop(k) - thresh) / max(ma(k), realmin) * wealth(k);
        s = (W - below - (wealth(k) - inside)) / W;
    end

mom.top1     = top_share(0.01);
mom.top5     = top_share(0.05);
mom.top10    = top_share(0.10);
mom.bottom50 = 1 - top_share(0.50);

% Gini from the Lorenz curve (trapezoid over the atomic distribution)
L = cum_wealth / W;
mom.gini = 1 - sum(ma .* (L + [0; L(1:end-1)]));

% Bequest flow and revenue
B = sum(g.aa .* sol.ops.death_rate .* m);
mom.bwr        = B / W;
mom.revenue    = sum(sol.ops.tax_fn(g.aa, params) .* sol.ops.death_rate .* m);
mom.entry_flow = sum(sol.ops.death_rate .* m);

% Profiles
mass_h             = accumarray(g.ih, m, [g.Nh, 1]);
mom.wealth_by_age  = accumarray(g.ih, g.aa .* m, [g.Nh, 1]) ./ mass_h;
mom.wealth_marginal = ma;

end
