function params = params_default()
% PARAMS_DEFAULT  Default parameter struct, Step 1 values (appendix section 8).
%
%   params = params_default()
%
% Single source of truth for all parameter values and the Step 1 <-> Step 2
% module toggles. Block modules (income_*.m, returns_*.m, ...) read everything
% from this struct; no parameter is hardcoded in module files.

%% Module toggles (operator_build.m dispatches on these)
params.income    = 'step1';  % 'step1': OU log income | 'step2': age profile + perm/trans
params.returns   = 'step1';  % 'step1': constant r0   | 'step2': heterogeneous returns
params.mortality = 'step1';  % 'step1': constant hazard | 'step2': Destatis spline
params.bequest   = 'step1';  % 'step1': warm glow | 'step2': De Nardi luxury form
params.tax       = 'step1';  % 'step1': flat above exemption | 'step2': full ErbStG

%% Preferences (Block F)
params.gamma = 2.0;          % CRRA risk aversion
params.rho   = 0.04;         % discount rate (annual)

%% Income process, Step 1 (Block A.1): dy = -theta_y (y - mu_y) dt + sigma_y dW
params.theta_y = 0.05;       % OU mean reversion of log income
params.sigma_y = 0.20;       % OU innovation std
params.mu_y    = 0.0;        % OU mean of log income
params.Ny      = 7;          % discrete income states
params.w       = 1.0;        % wage scale (normalisation)

%% Returns, Step 1 (Block B.1)
params.r0 = 0.03;            % safe real return

%% Mortality, Step 1 (Block C.1)
params.lambda_bar = 0.02;    % constant adult hazard (= 1/50)

%% Ages
params.h0    = 18;           % entry age
params.hR    = 65;           % retirement age
params.h_max = 99;           % terminal age (hard truncation)
params.hG    = 18;           % grant age (Block J: single lump sum at entry)

%% Bequest motive, Step 1 (Block D.1): W(a) = theta_b * a^(1-gamma)/(1-gamma)
params.theta_b = 1.0;        % warm-glow strength; PLACEHOLDER until calibrate_step1.m

%% Estate tax, Step 1 (Block E.1): T_e(b) = tau0 * max(b - F, 0)
params.tau0 = 0.20;          % flat rate (validation placeholder, not ErbStG)
params.F    = 4e5;           % exemption (EUR)

%% Capital grant (Block J)
params.G = 0;                % grant amount (EUR): 0 status quo, 2e4 Grunderbe

%% Grids (appendix section 3.1)
params.Na        = 300;      % wealth grid points
params.a_max     = 5e6;      % wealth upper bound (EUR)
params.a_curv    = 2;        % power-spacing curvature; > 1 refines near a = 0
params.Nh        = params.h_max - params.h0 + 1;   % = 82 age points, dh = 1
params.y_sd_span = 3;        % y grid spans mu_y +/- span * stationary OU sd

%% Numerics
params.Delta_hjb = 1000;     % implicit time step in HJB iteration
params.tol       = 1e-6;     % HJB convergence tolerance
params.maxit_hjb = 100;      % expected to converge within ~100 iterations

end
