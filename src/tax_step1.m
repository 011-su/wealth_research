function T = tax_step1(b, params)
% TAX_STEP1  Block E.1: flat estate-tax rate above a flat exemption.
%
%   T = tax_step1(b, params)
%
% Returns T, the estate tax due on bequest b (vectorised, same size as b):
%   T_e(b) = tau0 * max(b - F, 0),   b and F in model units.
%
% Unlike the other block modules this maps bequest values, not the state
% grid: kfe_solve.m's inheritance kernel calls it to compute heirs' initial
% wealth max(b - T_e(b), 0) + G, and the revenue integration calls it on
% the wealth grid. Step 1 rates are validation placeholders, not ErbStG.

T = params.tau0 * max(b - params.F, 0);

end
