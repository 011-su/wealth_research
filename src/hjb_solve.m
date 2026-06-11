function [V, c, adrift, A, n_iter] = hjb_solve(grids, ops, params)
% HJB_SOLVE  Implicit upwind finite-difference HJB solver (appendix 3.2).
%
%   [V, c, adrift, A, n_iter] = hjb_solve(grids, ops, params)
%
% Achdou et al. (2022) scheme, following Moll's huggett_diffusion_partialeq.m:
% upwind the wealth drift with the state-constraint boundary conditions
% V_a(0) = u'(c0) and V_a(a_max) = u'(c0), implicit update with time step
% Delta_hjb. The static generator blocks (income diffusion, aging, mortality)
% and the bequest RHS term come pre-assembled from operator_build.m; only the
% wealth-drift block A_a is rebuilt each iteration.
%
% Outputs: V, c, adrift (N x 1) value, consumption, optimal wealth drift;
%          A (N x N) converged generator (excl. the death->entry kernel);
%          n_iter iterations used.

N  = grids.N;
Na = grids.Na;
ga = params.gamma;

u     = @(x) x.^(1 - ga) / (1 - ga);
u1    = @(x) x.^(-ga);
u1inv = @(x) x.^(-1 / ga);

% Zero-saving resources and consumption (upwind fallback and boundaries)
res0 = ops.r .* grids.aa + ops.inc;
Va0  = u1(res0);

% Nonuniform forward/backward wealth steps per state (dummy at the edges,
% masked below)
da_f = grids.da(min(grids.ia, Na - 1));
da_b = grids.da(max(grids.ia - 1, 1));
up = grids.ia < Na;    % states with a forward neighbour
dn = grids.ia > 1;     % states with a backward neighbour
idx = (1:N)';

V = u(res0) / params.rho;   % initial guess: value of consuming res0 forever

for n = 1:params.maxit_hjb
    % One-sided differences with state-constraint boundary conditions
    Vaf = Va0;   % at a_max: forces sf = 0, no flux out of the grid
    Vab = Va0;   % at a = 0: V_a = u'(res0), borrowing constraint
    Vaf(up) = (V(idx(up) + 1) - V(up)) ./ da_f(up);
    Vab(dn) = (V(dn) - V(idx(dn) - 1)) ./ da_b(dn);
    Vaf = max(Vaf, 1e-12);   % guard against transient non-concavity
    Vab = max(Vab, 1e-12);

    cf = u1inv(Vaf);  sf = res0 - cf;
    cb = u1inv(Vab);  sb = res0 - cb;

    If = sf > 0;                % positive drift: forward difference
    Ib = sb < 0;                % negative drift: backward difference
    I0 = 1 - If - Ib;           % stationary: zero-saving consumption
    Va_up = Vaf .* If + Vab .* Ib + Va0 .* I0;

    c = u1inv(Va_up);

    % Wealth-drift block A_a (upwind, rows sum to zero)
    X = -min(sb, 0) ./ da_b;  X(~dn) = 0;
    Z =  max(sf, 0) ./ da_f;  Z(~up) = 0;
    Aa = sparse(idx, idx, -(X + Z), N, N) ...
       + sparse(idx(up), idx(up) + 1, Z(up), N, N) ...
       + sparse(idx(dn), idx(dn) - 1, X(dn), N, N);

    A = Aa + ops.Astatic;

    B   = (params.rho + 1 / params.Delta_hjb) * speye(N) - A;
    rhs = u(c) + V / params.Delta_hjb + ops.bequest_rhs;
    V_new = B \ rhs;

    dist = max(abs(V_new - V));
    V = V_new;
    if dist < params.tol
        n_iter = n;
        adrift = res0 - c;
        return
    end
end

warning('hjb_solve:noconvergence', ...
    'HJB did not converge in %d iterations (last dist %.2e)', params.maxit_hjb, dist);
n_iter = params.maxit_hjb;
adrift = res0 - c;

end
