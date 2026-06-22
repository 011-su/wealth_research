function out = solve_grant_flat_tax(params, G_eur, varargin)
% SOLVE_GRANT_FLAT_TAX  Grant funded by a flat surtax on ALL bequests (no exemption).
%
%   out = solve_grant_flat_tax(params, G_eur, 'Name', value, ...)
%
% Variant of Experiment 1 in which the grant is funded by a flat surtax
% tau_add on the WHOLE bequest (no exemption), levied on top of the unchanged
% status-quo estate tax (which stays as the held-constant leak). Roots tau_add
% so the surtax revenue funds the grant:
%
%     tau_add * E[b * death_rate]  =  G * N_entry.
%
% Because the surtax hits the whole bequest distribution (not just the thin
% slice above the €400k exemption), the base is far wider than in
% solve_revenue_balance, so much larger grants are feasible.
%
% As in solve_revenue_balance, the Step 1 HJB is invariant to the tax and the
% grant (gross-estate warm glow; grant enters only the entry kernel), so the
% HJB is solved ONCE and the root find re-solves only the KFE.
%
% Options:
%   'tau_add_max'  upper bound on the surtax (default 1 - tau0, so the top
%                  marginal rate on large estates stays <= 100%)
%   'verbose'      print progress (default true)
%
% Returns out: .feasible, .tau_add (NaN if infeasible), .sol/.mom at the
% solution, .sol_sq/.mom_sq (status quo, same grid/HJB), .surtax_rev, .bill,
% .entry_flow, .max_grant_eur (grant fundable at tau_add_max), .n_kfe.

ip = inputParser;
ip.addParameter('tau_add_max', 1 - params.tau0);
ip.addParameter('verbose', true);
ip.parse(varargin{:});
o = ip.Results;

if ~strcmp(params.bequest, 'step1')
    warning('solve_grant_flat_tax:invariance', ...
        'Solving one HJB for all tau_add assumes the Step 1 gross-estate warm glow.');
end

grids = grids_build(params);
ops   = operator_build(grids, params);
[V, c, adrift, A, n_iter] = hjb_solve(grids, ops, params);
if o.verbose, fprintf('HJB solved once (%d iters); sweeping tau_add (no-exemption surtax).\n', n_iter); end

n_kfe = 0; last = [];

    function s = eval_tg(tadd, Geur)
        pp = params; pp.G_eur = Geur; pp.tau_add = tadd; pp = params_derive(pp);
        tk = tic;
        [m, gd, info] = kfe_solve(grids, ops, pp, A);
        n_kfe = n_kfe + 1;
        s = struct('params', pp, 'grids', grids, 'ops', ops, 'A', A, ...
                   'V', V, 'c', c, 'adrift', adrift, 'm', m, 'g', gd, ...
                   'kfe', info, 'hjb_iter', n_iter);
        s.entry_flow = sum(ops.death_rate .* m);
        s.surtax_rev = tadd * sum(grids.aa .* ops.death_rate .* m);
        s.bill       = (Geur / params.eur_per_unit) * s.entry_flow;
        if o.verbose
            fprintf('  [kfe %d] tau_add=%.4f G=%.0fk -> surtax_rev %.5f vs bill %.5f (budget %+.2e, %.0f s)\n', ...
                n_kfe, tadd, Geur/1e3, s.surtax_rev, s.bill, s.surtax_rev - s.bill, toc(tk));
        end
    end
    % Root find evaluates at the target grant; bill uses G_eur.
    function f = resid(tadd), last = eval_tg(tadd, G_eur); f = last.surtax_rev - last.bill; end

% True status quo (G = 0, no surtax) for the comparison.
sol_sq = eval_tg(0, 0);

% Bracket the surtax: at tau_add = 0 with the grant present, surtax_rev = 0 < bill.
f_lo = resid(0);
f_hi = resid(o.tau_add_max);
out = struct(); out.tau_add_max = o.tau_add_max;
if f_hi < 0
    out.feasible = false; out.tau_add = NaN;
    out.max_grant_eur = last.surtax_rev / last.entry_flow * params.eur_per_unit;
else
    out.feasible = true;
    out.tau_add = fzero(@resid, [0, o.tau_add_max], optimset('TolX', 1e-5, 'Display', 'off'));
    resid(out.tau_add);
end

out.sol = last;   out.mom = moments(last, last.params);
out.sol_sq = sol_sq;  out.mom_sq = moments(sol_sq, sol_sq.params);
out.surtax_rev = last.surtax_rev;  out.bill = last.bill;  out.entry_flow = last.entry_flow;
out.n_kfe = n_kfe;

end
