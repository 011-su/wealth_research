function W = bequest_step1(grids, params)
% BEQUEST_STEP1  Block D.1: simple warm-glow bequest value.
%
%   W = bequest_step1(grids, params)
%
% Returns W (Na x 1): the warm-glow value of leaving bequest a,
%   W(a) = theta_b * (a + bequest_shift)^(1-gamma) / (1-gamma).
%
% The small shift keeps W(0) finite for gamma >= 1 (numerical
% regularisation, analogous to the De Nardi theta_2 shifter in Step 2).
% Enters the HJB as the inhomogeneous term lambda(h) * W(a).

W = params.theta_b * (grids.a + params.bequest_shift).^(1 - params.gamma) ...
    / (1 - params.gamma);

end
