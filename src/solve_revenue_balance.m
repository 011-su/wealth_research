function out = solve_revenue_balance(params, G_eur, varargin)
% SOLVE_REVENUE_BALANCE  Estate-tax rate funding a grant under incremental balance.
%
%   out = solve_revenue_balance(params, G_eur, 'Name', value, ...)
%
% Finds tau* such that the estate-tax revenue INCREMENT over the status quo
% funds the grant bill in the new stationary state (CLAUDE.md discrepancy 13):
%
%     Rev(tau*) = Rev_status_quo + G * N_entry.
%
% KEY SPEEDUP: in Step 1 the warm glow values the GROSS estate and the grant
% enters only the entry kernel, so the HJB generator A and the consumption
% policy are INVARIANT to tau and G. This routine therefore solves the HJB
% ONCE and the root find re-solves only the cheap stationary KFE per tau
% (one sparse solve each) via fzero. Cuts a ~6-solve outer loop from ~6 full
% HJB+KFE solves to 1 HJB + ~8 KFE solves.
%
% Name-value options:
%   'tau_max'  upper bound on the estate-tax rate (default 1.0 = 100% above
%              the exemption; the economic maximum)
%   'verbose'  print progress (default true)
%
% Returns struct out:
%   .feasible      true if the grant can be balanced with tau <= tau_max
%   .tau_star      balancing rate (NaN if infeasible)
%   .sol, .mom     full solution + moments at tau_star (or at tau_max if not)
%   .sol_sq, .mom_sq  status-quo solution + moments on the SAME grid/HJB
%                  (for self-consistent comparison figures, any resolution)
%   .rev_sq, .increment, .bill   revenue flow decomposition at the solution
%   .max_grant_eur approx. largest grant balanceable at tau_max (infeasible
%                  report; uses the revenue raised at tau_max)
%   .n_kfe         number of KFE solves used

ip = inputParser;
ip.addParameter('tau_max', 1.0);
ip.addParameter('verbose', true);
ip.parse(varargin{:});
o = ip.Results;

% Guard: the one-HJB-for-all-tau shortcut relies on the value function being
% invariant to the tax and grant. True for the Step 1 gross-estate warm glow;
% a Step 2 bequest motive over the NET estate would break it (re-solve per tau).
if ~strcmp(params.bequest, 'step1')
    warning('solve_revenue_balance:invariance', ...
        ['Reusing one HJB across tau assumes V is invariant to tau and G ' ...
         '(Step 1 gross-estate warm glow). For a net-of-tax bequest motive ' ...
         'this is INVALID -- re-solve the HJB per tau.']);
end

tau_sq = params.tau0;

% --- Solve the HJB ONCE (invariant to tau, G) ---
grids = grids_build(params);
ops   = operator_build(grids, params);
[V, c, adrift, A, n_iter] = hjb_solve(grids, ops, params);
if o.verbose
    fprintf('HJB solved once (%d iters); sweeping tau with KFE-only re-solves.\n', n_iter);
end

n_kfe = 0;
last  = [];

    function s = eval_tau(tau, Geur)
        pp = params; pp.tau0 = tau; pp.G_eur = Geur; pp = params_derive(pp);
        [m, gdens, info] = kfe_solve(grids, ops, pp, A);
        n_kfe = n_kfe + 1;
        s = struct('params', pp, 'grids', grids, 'ops', ops, 'A', A, ...
                   'V', V, 'c', c, 'adrift', adrift, 'm', m, 'g', gdens, ...
                   'kfe', info, 'hjb_iter', n_iter);
        s.entry_flow = sum(ops.death_rate .* m);
        s.revenue    = sum(ops.tax_fn(grids.aa, pp) .* ops.death_rate .* m);
    end

% Status-quo equilibrium (same HJB; status-quo rate, no grant). Kept for a
% grid-consistent comparison so the wrappers never depend on a saved
% status_quo.mat (which may be at a different resolution).
sol_sq = eval_tau(tau_sq, 0);
rev_sq = sol_sq.revenue;

G_model = G_eur / params.eur_per_unit;

    function f = resid(tau)
        last = eval_tau(tau, G_eur);
        f = last.revenue - rev_sq - G_model * last.entry_flow;
        if o.verbose
            fprintf('  tau = %.5f -> increment %.5f, bill %.5f, budget %+.2e (kfe %d)\n', ...
                tau, last.revenue - rev_sq, G_model * last.entry_flow, f, n_kfe);
        end
    end

f_lo = resid(tau_sq);
f_hi = resid(o.tau_max);

out = struct();
out.rev_sq  = rev_sq;
out.tau_max = o.tau_max;

if f_lo >= 0
    % Grant already funded at the status-quo rate (tiny grant)
    out.feasible = true; out.tau_star = tau_sq; resid(tau_sq);
elseif f_hi < 0
    % Infeasible: even tau_max cannot fund the grant
    out.feasible = false; out.tau_star = NaN;
    out.max_grant_eur = (last.revenue - rev_sq) / last.entry_flow * params.eur_per_unit;
else
    % Bracketed: hand to fzero (Brent). Evals are cheap KFE solves now.
    out.feasible = true;
    out.tau_star = fzero(@resid, [tau_sq, o.tau_max], ...
                         optimset('TolX', 1e-5, 'Display', 'off'));
    resid(out.tau_star);   % leave `last` at the solution
end

out.sol        = last;
out.mom        = moments(last, last.params);
out.sol_sq     = sol_sq;
out.mom_sq     = moments(sol_sq, sol_sq.params);
out.increment  = last.revenue - rev_sq;
out.bill       = G_model * last.entry_flow;
out.entry_flow = last.entry_flow;
out.n_kfe      = n_kfe;

if o.verbose
    fprintf('done: %d KFE solves, 1 HJB solve.\n', n_kfe);
end

end
