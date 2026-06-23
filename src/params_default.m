function params = params_default()
% PARAMS_DEFAULT  Default parameter struct, Step 1 values (appendix section 8).
%
%   params = params_default()
%
% Single source of truth for all parameter values and the Step 1 <-> Step 2
% module toggles. Block modules (income_*.m, returns_*.m, ...) read everything
% from this struct; no parameter is hardcoded in module files.
%
% Units: one model unit = eur_per_unit EUR (mean annual gross earnings), so
% w = 1 stays a normalisation while EUR-denominated policy parameters (a_max,
% F, G, ...) are converted below. Keeps c, V well scaled for the solvers.

%% Module toggles (operator_build.m dispatches on these)
params.income    = 'step1';  % 'step1': OU log income | 'step2': age profile + perm/trans
params.returns   = 'step1';  % 'step1': constant r0   | 'step2': heterogeneous returns
params.mortality = 'step1';  % 'step1': constant hazard | 'step2': Destatis spline
params.bequest   = 'step1';  % 'step1': warm glow | 'step2': De Nardi luxury form
params.tax       = 'step1';  % 'step1': flat above exemption | 'step2': full ErbStG

%% Units
params.eur_per_unit = 4.5e4; % EUR per model unit (approx. German mean annual
                             % gross earnings; Step 1 placeholder)

%% Preferences (Block F)
params.gamma = 2.0;          % CRRA risk aversion
params.rho   = 0.04;         % discount rate (annual)

%% Income process, Step 1 (Block A.1): dy = -theta_y (y - mu_y) dt + sigma_y dW
params.theta_y = 0.05;       % OU mean reversion of log income
params.sigma_y = 0.20;       % OU innovation std
params.mu_y    = 0.0;        % OU mean of log income
params.Ny      = 7;          % discrete income states
params.w       = 1.0;        % wage scale (normalisation)

% Retirement income: NOT in the appendix (budget has zero income for
% h >= hR, which makes c = 0 and u = -Inf at the constraint). Step 1
% placeholder at ~30% of mean earnings (Grundsicherung-like floor).
params.pension_eur = 1.35e4;

%% Returns, Step 1 (Block B.1)
params.r0 = 0.03;            % safe real return

%% Mortality, Step 1 (Block C.1)
params.lambda_bar = 0.02;    % constant adult hazard (= 1/50)

%% Ages
params.h0    = 18;           % entry age
params.hR    = 65;           % retirement age
params.h_max = 99;           % terminal age (hard truncation: aging out of
                             % h_max is treated as certain death with bequest)
params.hG    = 18;           % grant age (Block J: single lump sum at entry)

%% Bequest motive, Step 1 (Block D.1): W(a) = theta_b * (a + shift)^(1-gamma)/(1-gamma)
params.theta_b = 1.0;        % warm-glow strength; PLACEHOLDER until calibrate_step1.m
params.target_bwr = 0.024;   % calibration target: bequest flow / wealth per
                             % year. PLACEHOLDER: with the constant Step 1
                             % hazard the model caps bwr near lambda_bar
                             % (~0.020-0.026), so the German ~0.03 is only
                             % attainable with age-rising mortality
                             % (Step 2.C) -- see TODO.md
params.bequest_shift_eur = 2e4;  % small shift: W(0) finite for gamma >= 1
                             % (numerical regularisation, cf. De Nardi theta_2)

%% Estate tax, Step 1 (Block E.1): T_e(b) = tau0*max(b-F,0) + tau_add*b
params.tau0    = 0.20;       % flat rate above exemption (placeholder, not ErbStG)
params.F_eur   = 4e5;        % exemption (EUR)
params.tau_add = 0;          % flat grant-funding surtax on ALL bequests (no
                             % exemption); set by solve_grant_flat_tax.m, 0 = SQ

%% Capital grant (Block J)
params.G_eur = 0;            % grant (EUR): 0 status quo, 2e4 Grunderbe

%% Grids (appendix section 3.1)
params.Na        = 300;      % wealth grid points
params.a_max_eur = 5e6;      % wealth upper bound (EUR)
params.a_curv    = 2;        % power-spacing curvature; > 1 refines near a = 0
params.y_sd_span = 3;        % y grid spans mu_y +/- span * stationary OU sd

%% Numerics
params.Delta_hjb = 1000;     % implicit time step in HJB iteration
params.tol       = 1e-6;     % HJB convergence tolerance
params.maxit_hjb = 100;      % expected to converge within ~100 iterations
params.kfe_method = 'direct';% KFE solve: 'direct' (sparse LU) | 'iterative'
                             % (bicgstab + ilu0; bounded memory, no kernel fill-in)

%% Derived model-unit quantities (recompute via params_derive after
%% overriding any *_eur input)
params = params_derive(params);

end
