function time_kfe_methods(params)
% TIME_KFE_METHODS  Benchmark iterative (bicgstab+ilu0) vs direct (sparse LU)
% KFE solve, at the current grid, for G = 20k and G = 100k.
%
%   time_kfe_methods([params])
%
% Solves the HJB once (shared), then times both KFE solvers on the same
% (A' + R) system for each grant, and checks the two solutions agree.

if nargin < 1, params = params_default(); end

grids = grids_build(params);
ops   = operator_build(grids, params);
fprintf('Na=%d, N=%d. Solving HJB once (shared across solves)...\n', params.Na, grids.N);
[~, ~, ~, A] = hjb_solve(grids, ops, params);

fprintf('\n  G      | direct (sparse LU)     | iterative (gmres + ilu droptol)                 | speedup | max|dm|\n');
fprintf(  '---------+------------------------+-------------------------------------------------+---------+--------\n');
for Geur = [2e4, 1e5]
    p = params; p.G_eur = Geur; p = params_derive(p);

    p.kfe_method = 'direct';
    t = tic; [m_d, ~, info_d] = kfe_solve(grids, ops, p, A); td = toc(t);

    p.kfe_method = 'iterative';
    t = tic; [m_i, ~, info_i] = kfe_solve(grids, ops, p, A); ti = toc(t);

    err = max(abs(m_d - m_i));
    fprintf('  %4.0fk  | %6.1f s (res %.0e)   | %6.1f s %-6s %3.0f it res %.0e flag %d (ilu fill %.1fx) | %6.1fx | %.1e\n', ...
        Geur/1e3, td, info_d.residual, ti, info_i.kfe_solver, info_i.kfe_iter, ...
        info_i.residual, info_i.kfe_flag, info_i.ilu_fill, td/ti, err);
end
end
