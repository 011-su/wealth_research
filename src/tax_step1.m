function T = tax_step1(b, params)
% TAX_STEP1  Block E.1: flat estate-tax rate above a flat exemption.
%
%   T = tax_step1(b, params)
%
% Will return T, the estate tax due on bequest b (vectorised, same size
% as b):
%   T_e(b) = tau0 * max(b - F, 0).
%
% Unlike the other block modules this maps bequest values, not the state
% grid: kfe_solve.m's inheritance kernel calls it to compute heirs'
% initial wealth max(b - T_e(b), 0) + G, and equilibrium.m calls it for
% revenue integration. Step 1 rates are validation placeholders, not the
% ErbStG schedule (that is tax_step2.m).
%
% Stub (session 1) -- no implementation yet.

error('tax_step1: not implemented yet (scaffold stub)');

end
