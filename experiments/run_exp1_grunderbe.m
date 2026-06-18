function out = run_exp1_grunderbe(params)
% RUN_EXP1_GRUNDERBE  Experiment 1: Grunderbe, G = 20,000 EUR (outline §1).
%
%   out = run_exp1_grunderbe([params])
%
% Canonical Grunderbe: 20,000 EUR at entry (h = 18), funded by an estate-tax
% increment under incremental revenue balance (CLAUDE.md disc. 13). All logic
% (one HJB solve + fzero over KFE re-solves, status quo, reporting, figures)
% is in run_grunderbe.m / solve_revenue_balance.m. The 20,000 EUR here is the
% policy definition; nothing downstream hardcodes it. Larger-grant variant:
% run_exp1b_grunderbe_200k.m.

if nargin < 1, params = []; end
out = run_grunderbe(params, 2e4);

end
