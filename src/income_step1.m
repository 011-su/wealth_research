function [drift_y, var_y] = income_step1(grids, params)
% INCOME_STEP1  Block A.1: OU process on log income.
%
%   [drift_y, var_y] = income_step1(grids, params)
%
% Returns, on the y grid (Ny x 1):
%   drift_y  drift of y:     -theta_y * (y - mu_y)
%   var_y    variance of y:   sigma_y^2 (constant vector)
%
% The state is log income, so the OU coefficients enter directly (no Ito
% correction); operator_build.m assembles the upwinded diffusion generator
% L_y from these (same scheme as Moll's huggett_diffusion_partialeq.m).

drift_y = -params.theta_y * (grids.y - params.mu_y);
var_y   = params.sigma_y^2 * ones(grids.Ny, 1);

end
