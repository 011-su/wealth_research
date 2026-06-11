function [drift_y, var_y] = income_step1(grids, params)
% INCOME_STEP1  Block A.1: OU process on log income.
%
%   [drift_y, var_y] = income_step1(grids, params)
%
% Will return, on the y grid (Ny x 1):
%   drift_y  drift of y:     -theta_y * (y - mu_y)
%   var_y    variance of y:   sigma_y^2 (constant vector)
%
% operator_build.m uses these coefficients to assemble the upwinded
% diffusion generator L_y as a sparse Kronecker block in the y dimension
% (same scheme as Moll's huggett_diffusion_partialeq.m, but with the state
% kept in logs so the OU coefficients enter directly, no Ito correction).
%
% Stub (session 1) -- no implementation yet.

error('income_step1: not implemented yet (scaffold stub)');

end
