function r = returns_step1(grids, params)
% RETURNS_STEP1  Block B.1: single safe real rate.
%
%   r = returns_step1(grids, params)
%
% Will return r (N x 1): the return rate at every flattened state point.
% Step 1 specification: r(a) = r0 everywhere.
%
% operator_build.m uses r in the wealth drift
%   adot = r .* a + w * exp(y) .* (h < hR) - c.
%
% Stub (session 1) -- no implementation yet.

error('returns_step1: not implemented yet (scaffold stub)');

end
