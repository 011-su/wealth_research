# Implementation Appendix: continuous-time HA life-cycle model with bequests (MATLAB)

**Purpose.** Specify the mathematical model, the two-step implementation strategy, the numerical method, and the code structure for the preliminary paper. This document is the handoff specification for code generation.
**Companion document.** `preliminary_paper_outline.md` (paper structure and policy experiments).
**Implementation stack.** MATLAB R2021a+. Built on two existing libraries:

- **Moll reference codes** (<https://benjaminmoll.com/codes/>): the canonical implementation of the Achdou-Han-Lasry-Lions-Moll (2022) framework. Use as starting point for the HJB+KFE finite-difference machinery.
- **SparseEcon** (<https://github.com/schaab-lab/SparseEcon>): Schaab-Zhang toolbox providing cleaner abstractions for transition dynamics, parameter sweeps, and policy experiments, with enumerated use cases in `use_cases/`. Built on top of Moll's codebase.

The strategy is to inherit working HJB+KFE solvers from these libraries and add the project-specific structure (life-cycle aging, warm-glow bequests, the German estate-tax schedule, the inheritance kernel that links generations).

## 1. Model

### 1.1 State, controls, primitives

- State $x = (a, y, h)$: wealth $a \in [0, \bar a]$, log productivity $y \in [\underline y, \overline y]$, age $h \in [h_0, h_{\max}]$.
- Control: consumption $c \geq 0$.
- Primitives: discount rate $\rho$, risk aversion $\gamma$, wage scale $w$, retirement age $h_R$, terminal age $h_{\max}$, grant amount $G$, grant age $h_G$.
- Borrowing constraint: $a \geq 0$.

### 1.2 Dynamics

$$\dot a = r(a)\,a + w \cdot \exp(y) \cdot \mathbf{1}\{h < h_R\} + G \cdot \delta(h - h_G) - c,$$

$$dy = \mu_y(y, h)\, dt + \sigma_y(y, h)\, dW_t,$$

$$\dot h = 1.$$

At death (Poisson hazard $\lambda(h)$) the household leaves bequest $a$. The heir enters at $h = h_0$ with initial wealth $a_{\text{init}} = \max(a - T_e(a),\, 0) + G \cdot \mathbf{1}\{\text{grant active}\}$.

### 1.3 Preferences and bequest motive

$$u(c) = \frac{c^{1-\gamma}}{1-\gamma}, \qquad W(a) = \text{bequest value function (see Block D below)}.$$

### 1.4 HJB equation

$$\rho V(a,y,h) = \max_c \Big\{ u(c) + V_a \cdot \dot a + \mathcal{L}_y V + V_h \Big\} - \lambda(h)\big[V - W(a)\big],$$

with first-order condition $c = (V_a)^{-1/\gamma}$ and standard upwind treatment of $V_a$ at the borrowing constraint.

### 1.5 KFE equation

$$\partial_t g = -\partial_a(\dot a^\star\, g) + \mathcal{L}_y^\star g - \partial_h g - \lambda(h)\, g + B[g],$$

where $\dot a^\star$ is computed at the optimal policy and $B[\cdot]$ is the boundary operator injecting mass at $h = h_0$ to balance the mortality outflow. Stationary equilibrium: $\partial_t g = 0$.

## 2. Two-step implementation strategy

The model has six parameter blocks. Three are staged across Step 1 (working version, literature parameters; many lifted directly from Moll reference scripts) and Step 2 (realistic German calibration); three are single-step because there is no useful simpler alternative.

### Block A: Income process (*two-step*)

**Step 1.** Continuous-time OU on log income, no age dependence:

$$dy = -\theta_y (y - \mu_y)\, dt + \sigma_y\, dW_t.$$

Parameters from Drechsel-Grau et al. (2022) and Fuchs-Schündeln et al. (2010), at annual frequency:

$$\theta_y = 0.05, \quad \sigma_y = 0.20, \quad \mu_y = 0.$$

Discretise on $N_y = 7$ states. This matches the diffusion-income setup in Moll's `huggett_diffusion_partialeq.m` style, so the operator $\mathcal{L}_y$ is built using the same upwinded-diffusion scheme.

**Step 2.** Age-dependent productivity profile $\psi(h)$ plus persistent + transitory decomposition:

$$y_{i,h} = \psi(h) + \varepsilon^{\text{p}}_{i,h} + \varepsilon^{\text{t}}_{i,h},$$

with $\varepsilon^{\text{p}}$ AR(1) (same OU parameters as Step 1) and $\varepsilon^{\text{t}}$ i.i.d. Gaussian with $\sigma_{\text{t}} = 0.15$. The profile $\psi(h)$ is the German age-earnings profile from SOEP (hump-shaped, peak around age 50).

### Block B: Wealth returns (*two-step*)

**Step 1.** Single safe real rate:

$$r(a) \equiv r_0 = 0.03.$$

**Step 2.** Heterogeneous returns following Fagereng et al. (2020), adapted to PHF portfolio composition:

$$r(a, \omega) = r_0 + \beta \log\!\left(\frac{a}{\bar a_{\text{ref}}}\right) + \eta_\omega,$$

with persistent type $\omega$ drawn at $h = h_0$ from a discretised distribution (e.g. $N_\omega = 3$ types), $\beta \approx 0.005$ per log-EUR of wealth, $\eta_\omega$ cross-sectional dispersion with std 0.02 (annual). Reference wealth $\bar a_{\text{ref}}$ set to median household wealth from PHF.

### Block C: Mortality (*two-step*)

**Step 1.** Constant adult hazard:

$$\lambda(h) \equiv \bar\lambda = 1/50 = 0.02.$$

Expected adult life 50 years; $h_{\max}$ acts as a hard upper bound to truncate the tail.

**Step 2.** Age-specific hazard from Destatis *Sterbetafel* (latest available, sex-pooled):

$$\lambda(h) = \text{cubic spline fit to } \{q_h\}_{h = h_0}^{h_{\max}},$$

where $q_h$ is the one-year mortality rate at age $h$.

### Block D: Bequest motive (*two-step*)

**Step 1.** Simple warm-glow:

$$W(a) = \theta_b \cdot \frac{a^{1-\gamma}}{1 - \gamma}.$$

Single parameter $\theta_b$ calibrated to match the aggregate bequest-to-wealth ratio in the status-quo equilibrium.

**Step 2.** De Nardi (2004) luxury form:

$$W(a) = \theta_1 \cdot \frac{\left(1 + a / \theta_2\right)^{1-\gamma}}{1 - \gamma}.$$

Two parameters $(\theta_1, \theta_2)$ calibrated to jointly match the aggregate bequest-to-wealth ratio and the bequest-share-by-recipient-decile profile (Westermeier et al. 2016).

### Block E: Estate-tax schedule (*two-step*)

**Step 1.** Flat rate above flat exemption:

$$T_e(b) = \tau_0 \cdot \max(b - F, 0).$$

Placeholder values for code validation only: $\tau_0 = 0.20$, $F = 400{,}000$ EUR. These do not match the ErbStG schedule.

**Step 2.** Full *Erbschaftsteuergesetz* schedule:

- Three *Steuerklassen* (I: spouse/children/grandchildren; II: parents/siblings; III: all others) with class-specific progressive rate brackets and class-specific exemptions.
- *Verschonungsregeln* for business wealth: 85% reduction up to 26 M EUR (*Regelverschonung*); option for 100% (*Optionsverschonung*) under tighter conditions; sliding-scale phase-out above 26 M.
- Implemented as a piecewise-linear lookup on $b$ conditioned on (*Steuerklasse*, business-wealth share).

For the heir-class distribution, use Destatis fertility data combined with a simple rule (default: equal-share split among biological children of decedents with offspring; surviving spouse where applicable; *Steuerklasse* III for the childless-and-unmarried case).

### Single-step blocks

**Block F: Preferences.** CRRA throughout with $\gamma = 2$, $\rho = 0.04$.

**Block G: Numerical method.** Achdou et al. (2022) implicit upwind finite-difference HJB+KFE throughout, via Moll/SparseEcon. No two-step here.

**Block H: Borrowing constraint.** $a \geq 0$ throughout, enforced via boundary upwind condition $V_a(0, y, h) \geq u'(w \exp(y) \mathbf{1}\{h < h_R\})$.

**Block I: Equilibrium concept.** Partial equilibrium throughout, $r_0$ exogenous.

**Block J: Capital grant timing.** Lump-sum at $h_G$ throughout. The Frühstartrente accumulation (state contributions from age 6 to 18) is folded into a single equivalent lump sum at $h = 18$ given exogenous portfolio returns over the accumulation window.

## 3. Numerical method (MATLAB / Moll / SparseEcon)

### 3.1 Grid

- $a$: $N_a = 300$ points. Geometric spacing on $[0, \bar a]$ with $\bar a = 5 \times 10^6$ EUR. Refinement near $a = 0$ matters for the borrowing-constraint region.
- $y$: $N_y = 7$ (Step 1) or $N_y = 7 \times N_t$ (Step 2 with transitory shock superimposed).
- $h$: $N_h = h_{\max} - h_0 + 1$, uniform spacing $\Delta h = 1$ year. Default $h_0 = 18$, $h_{\max} = 99$, so $N_h = 82$.
- Total state dimension $N = N_a \cdot N_y \cdot N_h \approx 1.7 \times 10^5$ in Step 1.

Grid construction follows Moll's conventions: column-major flattening of the $(a, y, h)$ tensor, sparse Kronecker-product operators built per dimension and combined additively.

### 3.2 HJB solver (use existing infrastructure)

Do not write the HJB solver from scratch. Adapt the closest Moll reference script (a stationary-equilibrium script with diffusion income; check `huggett_diffusion_partialeq.m` and the life-cycle scripts on Moll's codes page for the closest starting point). The extensions to make:

1. Add the age dimension $h$ as a third state. Deterministic forward upwind: $\partial_h V$ becomes a sub-diagonal block of the operator.
2. Add the death-with-bequest term: a diagonal contribution $-\lambda(h)$ on the LHS, plus an inhomogeneous term $\lambda(h) W(a)$ on the RHS.
3. Replace the income process / returns / bequest function with module functions (see code organisation below) so that Step 1 ↔ Step 2 swaps are localised.

The implicit upwind iteration follows the Achdou et al. scheme (cf. Moll's numerical appendix):

```matlab
% Pseudocode skeleton (adapted from huggett_diffusion_partialeq.m structure)
V = V0;                                 % initial guess
for n = 1:max_iter
    [Va_fwd, Va_bwd] = compute_Va(V, a_grid);  % forward / backward differences
    c_fwd = Va_fwd.^(-1/gamma);
    c_bwd = Va_bwd.^(-1/gamma);
    [c_star, Va_up] = upwind(c_fwd, c_bwd, a_grid, y_grid, h_grid, ...);
    A = build_operator(c_star, params);   % a-drift + L_y + h-aging + mortality
    LHS = (rho + 1/Delta) * speye(N) - A;
    RHS = u(c_star) + V/Delta + lambda_h .* W_of_a;
    V_new = LHS \ RHS;                    % sparse linear solve
    if max(abs(V_new - V), [], 'all') < tol; break; end
    V = V_new;
end
```

### 3.3 KFE solver

For stationary equilibrium, use the same operator $A$ transposed plus the boundary operator that injects mass at $h = h_0$. SparseEcon's stationary-equilibrium routines provide a clean wrapper; otherwise adapt directly from Moll's KFE code (the KFE in his stationary scripts is solved jointly with the HJB).

```matlab
% Pseudocode skeleton
AT = A';                                  % transpose of converged HJB operator
% Inject mass at h = h_0 to balance deaths: for each (a, y, h) cell with h > h_0,
% mass lambda(h)*g(a,y,h) leaves via death. Add inflow at h = h_0 with wealth
% distributed via inheritance kernel (post-tax bequest + grant).
AT_mod = inject_inheritance_kernel(AT, h0_idx, tax_fn, grant);
g = solve_stationary_KFE(AT_mod);        % null-space-of-AT_mod approach
g = g / sum(g(:) .* da .* dy .* dh);     % normalise mass to 1
```

For transition dynamics, use SparseEcon's transition-dynamics routine (`use_cases/05_...` or equivalent). The pattern is backward HJB sweep, forward KFE sweep, possibly iterated for general equilibrium (not needed here since we are in partial equilibrium).

### 3.4 Revenue balance (outer loop)

For Experiment 1 (and 2): bisect on the estate-tax schedule parameter $\tau$ to satisfy
$\int T_e \lambda g\, dx = G \cdot N_{\text{entry}}$ at the new stationary $g$. Inner loop = full HJB+KFE solve.

```matlab
% Outer loop: bisection on tau
tau_lo = 0; tau_hi = 0.6;
for outer = 1:50
    tau_mid = 0.5 * (tau_lo + tau_hi);
    params.tau = tau_mid;
    [V, g] = solve_stationary(params);
    rev = compute_revenue(g, params);
    budget = rev - G * Nentry;
    if abs(budget) < 1e-3 * G * Nentry; break; end
    if budget < 0; tau_lo = tau_mid; else; tau_hi = tau_mid; end
end
```

### 3.5 Transition dynamics

For transition from pre-reform $g_0$ to post-reform stationary $g^\star$:

1. Discretise time on $[0, T]$ with $T$ large enough (e.g. 200 years for steady-state-to-steady-state convergence with the warm-glow bequest motive).
2. Backward sweep: solve HJB from $t = T$ (terminal $V_T = V^\star$) back to $t = 0$, treating policy parameters as fixed at the post-reform values.
3. Forward sweep: solve KFE from $t = 0$ (initial $g_0$) to $t = T$ using the policies from step 2.
4. In partial equilibrium with no time-varying prices, no iteration is needed: a single backward + forward pass suffices.

SparseEcon's transition-dynamics template handles the bookkeeping. The key project-specific element is making sure the boundary operator (deaths recycled to entries) is applied correctly at every time step in the forward sweep.

## 4. Calibration procedure

### 4.1 Step 1

All parameters set directly to literature values listed above. Only $\theta_b$ in the warm-glow Step 1 specification is calibrated, via simple bisection to match the target bequest-to-wealth ratio in the Step 1 stationary equilibrium. Total Step 1 calibration runtime: a few HJB+KFE solves.

### 4.2 Step 2

Parameters that are still directly observable (mortality, age profile, return mean, return dispersion) are set externally. Parameters identified indirectly (the bequest-motive parameters $\theta_1, \theta_2$ in Step 2) are calibrated by indirect inference:

1. Compute target moments $\mathbf{m}^{\text{data}}$ from German data sources.
2. Define loss $\mathcal{L}(\theta) = \|\mathbf{m}^{\text{model}}(\theta) - \mathbf{m}^{\text{data}}\|^2_W$ with weighting matrix $W$.
3. Minimise via Nelder-Mead (`fminsearch`) or differential evolution; no analytic gradient. Budget $\approx 100$ HJB+KFE solves.

Target moments to be matched in Step 2 calibration: aggregate bequest-to-wealth ratio; bequest-share-by-recipient-decile; wealth share held by households over 65.

## 5. Validation checkpoints

### 5.1 Step 1 validation

1. HJB converges within $\sim 100$ iterations to tolerance $10^{-6}$.
2. KFE solution has $\int g\, dx = 1$ to machine precision.
3. **Aiyagari limit.** With $h$ collapsed to a single age and infinite horizon (no mortality, no aging), the model reduces to Aiyagari (1994) with diffusion income. Compare to Moll's stationary Aiyagari diffusion script directly — output should match within numerical tolerance. This is essentially a regression test against the upstream code.
4. **Life-cycle no-mortality limit.** Set $\lambda \to 0$ and the warm-glow weight $\theta_b \to 0$. The model collapses to a deterministic life-cycle problem with a hard terminal age. Compare $V$ and consumption policy to an analytic Ramsey-style solution where one is available, or to a separately implemented dynamic-programming benchmark.
5. **No-income-heterogeneity limit.** Set $\sigma_y = 0$. The model collapses to a deterministic life-cycle problem. Sanity check the consumption profile.
6. Mass-conservation check on the boundary operator: total inflow at $h = h_0$ from the inheritance kernel equals total mortality outflow.

### 5.2 Step 2 validation

1. Step 1 result reproducible by setting Step 2 toggles to Step 1 values (the module switches in `params` and the function pointers in the block files should make this trivial).
2. All Step 1 limits still hold.
3. Calibration converges to target moments within tolerance $\sim 5\%$.
4. ABS top-share series (out-of-calibration) lies within 2–3 percentage points of model prediction in the status-quo equilibrium.

### 5.3 Per-experiment validation

1. Revenue balance holds to within $0.1\%$ in the new stationary state.
2. Transition path connects $g_0$ to $g^\star$ monotonically (or, if non-monotone, the non-monotonicity is interpretable).
3. Total wealth changes by an amount consistent with the cumulative transfers in the transition.

## 6. Code organisation

MATLAB project layout. SparseEcon-style convention: scripts at the top level, model-specific helpers in a `src/` directory, parameters as a `params` struct passed everywhere.

```
preliminary_paper/
  README.md
  startup.m              % adds paths, sets formatting; called at session start
  src/
    params_default.m     % returns the default params struct
    grids_build.m        % constructs a / y / h grids and sparse operators
    income_step1.m       % Block A.1: OU income process generator
    income_step2.m       % Block A.2: age profile + perm/trans
    returns_step1.m      % Block B.1: constant safe rate
    returns_step2.m      % Block B.2: Fagereng-style heterogeneous returns
    mortality_step1.m    % Block C.1: constant hazard
    mortality_step2.m    % Block C.2: Destatis spline
    bequest_step1.m      % Block D.1: simple warm-glow value W(a)
    bequest_step2.m      % Block D.2: De Nardi luxury form
    tax_step1.m          % Block E.1: flat-rate-above-exemption
    tax_step2.m          % Block E.2: full ErbStG schedule
    operator_build.m     % builds the sparse HJB operator A from the block functions
    hjb_solve.m          % implicit upwind HJB solver
    kfe_solve.m          % stationary KFE solver with inheritance kernel
    transition_solve.m   % backward HJB + forward KFE for transition dynamics
    equilibrium.m        % wrapper: returns (V, g) for given params; handles bisection on tau
    moments.m            % computes top shares, Gini, bequest-to-wealth ratio, ...
    calibrate_step1.m    % bisects on theta_b to match BWR
    calibrate_step2.m    % indirect inference on (theta_1, theta_2)
  experiments/
    run_status_quo.m
    run_exp1_grunderbe.m
    run_exp2_tax_reform.m
    run_exp3_universal_transfer.m
    run_exp4_verschonung_closure.m
    run_all.m            % calls all five; saves results to results/
  tests/
    test_aiyagari_limit.m
    test_lifecycle_limit.m
    test_mass_conservation.m
    test_step1_reproducible_from_step2.m
  postprocess/
    make_table1.m        % calibration table
    make_table2.m        % summary policy results
    make_figures.m       % Figures 1-4 for the paper
  data/
    destatis_sterbetafel.csv
    soep_age_profile.csv
    targets_step2.mat    % calibration moments and weights
  external/
    Moll-codes/          % git submodule or local copy of relevant Moll scripts
    SparseEcon/          % git submodule of schaab-lab/SparseEcon
  results/
    step1/
    step2/
```

The `params` struct carries Step 1 ↔ Step 2 toggles: e.g. `params.income = 'step1'` selects `income_step1.m`. A small dispatcher in `operator_build.m` chooses the right module functions based on the toggles.

## 7. Development order

This order leverages existing Moll/SparseEcon infrastructure aggressively.

1. **Bootstrap.** Clone SparseEcon and Moll codes into `external/`. Run SparseEcon's `use_cases/01_huggett` and verify the toolchain works on the local MATLAB install. Identify the closest reference script (likely SparseEcon's Huggett-diffusion variant, or a life-cycle template if available; otherwise Moll's Huggett-diffusion partial-equilibrium script).
2. **Scaffold.** Create `params_default.m`, `grids_build.m`, and stub files for all Step 1 blocks. Run the existing reference script unchanged through your new scaffold to confirm nothing is broken.
3. **Add the age dimension** to the operator builder. Run the no-mortality, no-bequest limit and check that the wealth distribution by age looks reasonable (mean wealth rising with age then falling near retirement).
4. **Add mortality and bequests.** Implement `mortality_step1.m`, `bequest_step1.m`, and the inheritance-kernel boundary operator in `kfe_solve.m`. Validate mass conservation.
5. **Status-quo equilibrium (Step 1).** Implement `tax_step1.m` with placeholder rates, run `equilibrium.m`. Plot the wealth density. First sanity-check result.
6. **Bisection on $\theta_b$** to match a target bequest-to-wealth ratio. Implement `calibrate_step1.m`.
7. **Revenue-balance outer loop.** Implement Experiment 1 (`run_exp1_grunderbe.m`). Bisect on $\tau$ for revenue balance.
8. **Experiments 2-4.** Mostly variations on Experiment 1.
9. **Transition dynamics.** Implement `transition_solve.m`, leveraging SparseEcon's transition routines. Re-run all experiments with full transition paths.
10. **Step 2.C (mortality)** — smallest change. Add `mortality_step2.m` and parameterise Destatis spline.
11. **Step 2.A (income process)** — add age profile, perm/trans decomposition.
12. **Step 2.B (returns)** — add heterogeneous return types.
13. **Step 2.D (bequest motive)** — add De Nardi form. Re-calibrate via indirect inference (`calibrate_step2.m`).
14. **Step 2.E (full ErbStG)** — implement Steuerklassen and Verschonungsregeln.
15. **Re-run all experiments** with Step 2 calibration. Generate paper Tables 1-2 and Figures 1-4 via `postprocess/`.
16. **Sensitivity runs** for the discussion section.

After Step 7 you have a publishable Step 1 set of preliminary results — useful for an early presentation or seminar slide deck even before Step 2 is done.

## 8. Reference parameter table (Step 1)

| Symbol | Value | Description |
|---|---|---|
| `gamma` | 2.0 | CRRA risk aversion |
| `rho` | 0.04 | Discount rate (annual) |
| `theta_y` | 0.05 | OU mean-reversion of log income |
| `sigma_y` | 0.20 | OU innovation std |
| `mu_y` | 0.0 | OU mean of log income |
| `Ny` | 7 | Discrete income states |
| `r0` | 0.03 | Safe real return |
| `lambda_bar` | 0.02 | Constant adult mortality hazard |
| `h0` | 18 | Entry age |
| `hR` | 65 | Retirement age |
| `h_max` | 99 | Terminal age |
| `w` | 1.0 | Wage scale (normalisation) |
| `theta_b` | calibrated | Warm-glow strength (Block D.1) |
| `tau0` | 0.20 | Step 1 flat estate-tax rate |
| `F` | 4e5 | Step 1 estate-tax exemption (EUR) |
| `G` | 0 / 2e4 | Capital grant (status quo / Grunderbe) |
| `Na` | 300 | Grid points in wealth |
| `a_max` | 5e6 | Wealth grid upper bound (EUR) |
| `Nh` | 82 | Grid points in age |
| `Delta_hjb` | 1000 | Implicit time step |
| `tol` | 1e-6 | HJB convergence tolerance |

## 9. Key external references for Claude Code

When you need to look up implementation details:

- **Moll's online numerical appendix:** <https://benjaminmoll.com/wp-content/uploads/2020/02/HACT_Numerical_Appendix.pdf> — definitive reference on the upwind scheme, sparse operator construction, borrowing-constraint handling.
- **Moll codes page:** <https://benjaminmoll.com/codes/> — list of scripts by topic. The "Stationary Eq. with Diffusion" and the life-cycle entries are the closest starting points.
- **SparseEcon README and use_cases:** <https://github.com/schaab-lab/SparseEcon> — abstraction layer and reusable templates.
- **Achdou et al. (2022), *Review of Economic Studies* 89(1).** The published paper, definitive for the mathematical setup.
