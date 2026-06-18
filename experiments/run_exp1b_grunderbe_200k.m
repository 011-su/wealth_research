function out = run_exp1b_grunderbe_200k(params)
% RUN_EXP1B_GRUNDERBE_200K  Large-grant variant of Experiment 1: G = 200,000 EUR.
%
%   out = run_exp1b_grunderbe_200k([params])
%
% Tenfold grant relative to exp1, to stress the estate-tax base: a grant to
% every entrant is funded only from estates above the exemption, so a large
% grant may be INFEASIBLE under revenue balance. run_grunderbe.m detects and
% reports that case (with the largest balanceable grant). The 200,000 EUR here
% is the policy definition; nothing downstream hardcodes it.

if nargin < 1, params = []; end
out = run_grunderbe(params, 2e5);

end
