function lambda = mortality_step1(grids, params)
% MORTALITY_STEP1  Block C.1: constant adult mortality hazard.
%
%   lambda = mortality_step1(grids, params)
%
% Will return lambda (Nh x 1): the mortality hazard at each age.
% Step 1 specification: lambda(h) = lambda_bar for all h.
%
% operator_build.m places -lambda(h) on the operator diagonal (death
% outflow) and pairs it with the inhomogeneous bequest term
% lambda(h) * W(a) on the HJB right-hand side; kfe_solve.m uses the same
% lambda for the mortality outflow balanced by the entry inflow at h = h0.
%
% Stub (session 1) -- no implementation yet.

error('mortality_step1: not implemented yet (scaffold stub)');

end
