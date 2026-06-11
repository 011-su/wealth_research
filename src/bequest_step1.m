function W = bequest_step1(grids, params)
% BEQUEST_STEP1  Block D.1: simple warm-glow bequest value.
%
%   W = bequest_step1(grids, params)
%
% Will return W (Na x 1): the warm-glow value of leaving bequest a,
%   W(a) = theta_b * a^(1-gamma) / (1-gamma).
%
% Enters the HJB as the inhomogeneous term lambda(h) * W(a) on the
% right-hand side (assembled in operator_build.m / hjb_solve.m).
% Note for the implementation session: W(0) = -Inf for gamma >= 1; needs
% the usual small-offset treatment at a = 0 (cf. v_terminal in Moll's
% lifecycle.m).
%
% Stub (session 1) -- no implementation yet.

error('bequest_step1: not implemented yet (scaffold stub)');

end
