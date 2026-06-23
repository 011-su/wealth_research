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
% The linear solve uses params.kfe_method: 'direct' (sparse LU, default) or
% 'iterative' (bicgstab + ilu(0) preconditioner). The iterative path avoids
% the LU fill-in that the inheritance kernel R induces, keeping memory bounded
% and (at fine grids / large grants) running far faster.
%
% Inputs:  A from hjb_solve (generator excluding the kernel).
% Outputs: m (N x 1) masses summing to 1; g (N x 1) density m ./ wx;
%          info: kernel matrix R, residual, column-sum check, entry flow,
%                and bicg_iter/bicg_flag/bicg_relres for the iterative path.

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

% Fix one row to pin the scale (Moll); M m = rhs is then nonsingular.
i_fix = 1;
rhs = zeros(N, 1); rhs(i_fix) = 0.1;
M = M0;
M(i_fix, :) = sparse(1, i_fix, 1, 1, N);

method = 'direct';
if isfield(params, 'kfe_method') && ~isempty(params.kfe_method)
    method = params.kfe_method;
end

info.kfe_solver = method; info.kfe_iter = NaN; info.kfe_flag = 0;
info.kfe_relres = NaN; info.ilu_fill = NaN;
switch method
    case 'direct'
        % Sparse LU (mldivide). Exact, but the inheritance kernel R couples
        % every dying state to the h0 block, so the fill-in -- and hence time
        % and memory -- blow up at fine grids / large grants.
        m = M \ rhs;
    case 'iterative'
        % Krylov solve with an incomplete-LU preconditioner: only matvecs +
        % triangular solves, so memory stays ~nnz of the (incomplete) factors,
        % far below the full LU fill-in the kernel R induces. ilu(0) is too
        % weak here (bicgstab breaks down), so use a small-drop-tolerance ilu
        % -- much stronger, still a small fraction of a full LU.
        [L, U] = ilu(M, struct('type', 'crout', 'droptol', 1e-4));
        info.ilu_fill = (nnz(L) + nnz(U)) / nnz(M);     % vs full LU this stays small
        % The BiCGSTAB family (bicgstab, bicgstabl) BREAKS DOWN here (flag 4)
        % on this nonsymmetric fix-one-row generator, regardless of
        % preconditioner -- verified at Na=100. Restarted gmres is robust and
        % converges in ~12 iterations; we try bicgstabl first (cheap, fails
        % fast) only to record the breakdown, then fall through to gmres.
        [m, flag, relres, iter] = bicgstabl(M, rhs, 1e-10, 500, L, U);
        info.kfe_solver = 'bicgstabl';
        if flag ~= 0
            [m, flag, relres, gi] = gmres(M, rhs, 100, 1e-10, 10, L, U);
            iter = (gi(1) - 1) * 100 + gi(2); info.kfe_solver = 'gmres';
        end
        info.kfe_iter = iter; info.kfe_flag = flag; info.kfe_relres = relres;
        % flag 3 = stagnation; benign once the residual is already ~1e-16.
        if flag ~= 0 && relres > 1e-6
            warning('kfe_solve:iter', '%s not converged: flag=%d relres=%.2e', ...
                info.kfe_solver, flag, relres);
        end
    otherwise
        error('kfe_solve:method', 'unknown params.kfe_method "%s"', method);
end
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
