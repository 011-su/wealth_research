function ops = operator_build(grids, params)
% OPERATOR_BUILD  Assemble the static parts of the HJB generator (appendix 3.1-3.2).
%
%   ops = operator_build(grids, params)
%
% Dispatches on the Step 1/Step 2 toggles in params and returns:
%   .Ay        (N x N)  income generator L_y: upwinded OU diffusion on the y
%                       grid, reflecting at both ends, Kronecker-lifted
%   .Ah        (N x N)  aging operator: deterministic unit drift, forward
%                       upwind; at h_max the outflow leaves the grid (aging
%                       out = certain death with bequest)
%   .Astatic   (N x N)  Ay + Ah - diag(lambda): everything in the generator
%                       that does not depend on the consumption policy
%   .r, .inc   (N x 1)  return rate; non-asset income w*exp(y) while working,
%                       pension after hR
%   .lambda    (N x 1)  mortality hazard by state
%   .exit      (N x 1)  aging-out rate (1/dh at h = h_max, else 0)
%   .death_rate(N x 1)  lambda + exit: total bequest-leaving outflow rate
%   .W         (N x 1)  warm-glow bequest value
%   .bequest_rhs (N x 1) death_rate .* W: inhomogeneous HJB RHS term
%   .p_y       (Ny x 1) stationary distribution of the discretised income
%                       process (entry distribution of heirs' y)
%   .tax_fn    handle   estate-tax function T(b, params)
%
% The policy-dependent wealth-drift block A_a is assembled per iteration in
% hjb_solve.m (and per time step in the future transition_solve.m); only the
% blocks above are static, so transitions can reuse this struct unchanged.

income_fn    = str2func(['income_'    params.income]);
returns_fn   = str2func(['returns_'   params.returns]);
mortality_fn = str2func(['mortality_' params.mortality]);
bequest_fn   = str2func(['bequest_'   params.bequest]);
ops.tax_fn   = str2func(['tax_'       params.tax]);

N  = grids.N;
Na = grids.Na;
Ny = grids.Ny;
Nh = grids.Nh;

% Returns and non-asset income
ops.r = returns_fn(grids, params);
working = grids.hh < params.hR;
ops.inc = params.w * exp(grids.yy) .* working + params.pension .* ~working;

% Mortality and bequest value, lifted to the flattened grid
lambda_h   = mortality_fn(grids, params);
ops.lambda = lambda_h(grids.ih);
W_a        = bequest_fn(grids, params);
ops.W      = W_a(grids.ia);

% Income generator on the y grid (Moll upwind scheme, reflecting barriers)
[drift_y, var_y] = income_fn(grids, params);
dy   = grids.dy;
chi  = -min(drift_y, 0) / dy + var_y / (2 * dy^2);   % rate down (y -> y-)
zeta =  max(drift_y, 0) / dy + var_y / (2 * dy^2);   % rate up   (y -> y+)
diag_y = -(chi + zeta);
diag_y(1)   = diag_y(1)   + chi(1);    % reflect blocked flux at the edges
diag_y(end) = diag_y(end) + zeta(end);
Ay1 = sparse([1:Ny, 1:Ny-1, 2:Ny], [1:Ny, 2:Ny, 1:Ny-1], ...
             [diag_y; zeta(1:Ny-1); chi(2:Ny)], Ny, Ny);
ops.Ay = kron(speye(Nh), kron(Ay1, speye(Na)));

% Stationary distribution of the discretised income process (heirs' entry y)
ops.p_y = [Ay1'; ones(1, Ny)] \ [zeros(Ny, 1); 1];

% Aging operator: unit drift, forward upwind; pure outflow at h_max
if Nh > 1
    dh  = grids.dh;
    Ah1 = spdiags([-ones(Nh, 1), ones(Nh, 1)] / dh, [0, 1], Nh, Nh);
    ops.Ah   = kron(Ah1, speye(Na * Ny));
    ops.exit = (grids.ih == Nh) / dh;
else
    ops.Ah   = sparse(N, N);   % age collapsed: infinite-horizon limit
    ops.exit = zeros(N, 1);
end

ops.Astatic     = ops.Ay + ops.Ah - spdiags(ops.lambda, 0, N, N);
ops.death_rate  = ops.lambda + ops.exit;
ops.bequest_rhs = ops.death_rate .* ops.W;

end
