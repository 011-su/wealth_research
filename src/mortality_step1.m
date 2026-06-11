function lambda = mortality_step1(grids, params)
% MORTALITY_STEP1  Block C.1: constant adult mortality hazard.
%
%   lambda = mortality_step1(grids, params)
%
% Returns lambda (Nh x 1): the mortality hazard at each age.
% Step 1 specification: lambda(h) = lambda_bar for all h. operator_build.m
% puts -lambda on the generator diagonal and pairs it with the bequest term
% lambda * W(a) on the HJB right-hand side; the same rates feed the
% death -> entry inheritance kernel in kfe_solve.m.

lambda = params.lambda_bar * ones(grids.Nh, 1);

end
