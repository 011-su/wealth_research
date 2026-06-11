function [m, g, info] = kfe_solve(grids, ops, params, A)
% KFE_SOLVE  Stationary KFE with the death -> entry inheritance kernel (appendix 3.3).
%
%   [m, g, info] = kfe_solve(grids, ops, params, A)
%
% Works in point masses m (m = g .* wx), for which the upwind generator is
% an exact transition-rate matrix: stationarity is (A' + R) m = 0 with R the
% inheritance kernel. R recycles every bequest-leaving flow -- mortality
% lambda(h) plus aging out at h_max -- into entry at h = h0: heirs start
% with wealth a' = max(b - T_e(b), 0) + G (linearly interpolated onto the a
% grid) and an income state drawn from the stationary distribution p_y.
% Solved by the fix-one-row trick, then normalised to unit mass.
%
% Inputs:  A from hjb_solve (generator excluding the kernel).
% Outputs: m (N x 1) masses summing to 1; g (N x 1) density m ./ wx;
%          info: kernel matrix R, residual, column-sum check, entry flow.

N = grids.N;

if any(ops.death_rate > 0)
    if params.hG ~= params.h0
        error('kfe_solve:grant_timing', ...
            'Only hG = h0 implemented (grant folded into the entry kernel).');
    end
    % Heir wealth as a function of the estate, computed once on the a grid
    b    = grids.a;
    anet = max(b - ops.tax_fn(b, params), 0) + params.G;
    anet = min(max(anet, grids.a(1)), grids.a(end));
    [ilo, wlo] = interp_weights(grids.a, anet);

    % Kernel: source (a, y, h) -> targets (a'(a), y', h0), y' ~ p_y
    tlo = ilo(grids.ia);          % N x 1 target a-index (lower neighbour)
    w1  = wlo(grids.ia);
    dr  = ops.death_rate;
    ii = zeros(2 * N * grids.Ny, 1); jj = ii; vv = ii;
    k = 0;
    for jy = 1:grids.Ny
        off = (jy - 1) * grids.Na;     % h0 block sits at offset 0
        ii(k+1 : k+2*N) = [tlo + off; tlo + 1 + off];
        jj(k+1 : k+2*N) = [(1:N)'; (1:N)'];
        vv(k+1 : k+2*N) = [dr .* w1; dr .* (1 - w1)] * ops.p_y(jy);
        k = k + 2 * N;
    end
    R = sparse(ii, jj, vv, N, N);
else
    R = sparse(N, N);   % no deaths (collapsed/limit configurations)
end

M0 = A' + R;

% Fix one row to pin the scale (Moll), then renormalise
i_fix = 1;
rhs = zeros(N, 1); rhs(i_fix) = 0.1;
M = M0;
M(i_fix, :) = sparse(1, i_fix, 1, 1, N);
m = M \ rhs;
m = m / sum(m);

g = m ./ grids.wx;

info.R          = R;
info.residual   = full(max(abs(M0 * m)));
info.colsum_max = full(max(abs(sum(M0, 1))));   % exact conservation: ~ 0
info.entry_flow = sum(ops.death_rate .* m); % total recycled mass per unit time

end

function [ilo, wlo] = interp_weights(grid, x)
% Linear interpolation weights of points x onto grid: x ~ wlo * grid(ilo)
% + (1-wlo) * grid(ilo+1), with ilo in [1, numel(grid)-1].
ilo = discretize(x, grid);
ilo = min(max(ilo, 1), numel(grid) - 1);
wlo = (grid(ilo + 1) - x) ./ (grid(ilo + 1) - grid(ilo));
wlo = min(max(wlo, 0), 1);
end
