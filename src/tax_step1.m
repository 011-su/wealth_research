function T = tax_step1(b, params)
% TAX_STEP1  Block E.1: flat estate tax above a flat exemption, plus an
% optional flat surtax on the WHOLE bequest (no exemption).
%
%   T = tax_step1(b, params)
%
%   T_e(b) = tau0 * max(b - F, 0)     [status-quo schedule, exemption F]
%          + tau_add * b               [grant-funding surtax, NO exemption]
%
% tau_add (params.tau_add, default 0) is the additional flat rate used by the
% no-exemption variant of Experiment 1 (solve_grant_flat_tax.m): the grant is
% funded by a surtax that hits every bequest from the first euro, while the
% status-quo schedule (with its exemption) is held constant as the leak. With
% tau_add = 0 this reduces to the plain status-quo schedule. b, F in model units.

tau_add = 0;
if isfield(params, 'tau_add'), tau_add = params.tau_add; end

T = params.tau0 * max(b - params.F, 0) + tau_add * b;

end
