# CLAUDE.md — WP0 preliminary paper implementation

Continuous-time HA life-cycle model with bequests, German estate-tax /
capital-grant experiments. Specification: `docs/implementation_appendix.md` +
`docs/preliminary_paper_outline.md` (the appendix is the contract). MATLAB,
building on Moll reference codes and SparseEcon (both in `external/`,
read-only, git-ignored).

## State after session 1 (2026-06-11)

- **Toolchain verified.** MATLAB R2026a Update 2 (Apple Silicon,
  `/Applications/MATLAB_R2026a.app`), headless `-batch` mode works, sparse
  matrices fine. Toolboxes: Optimization, Statistics/ML, Econometrics,
  Parallel, Signal. **Note: trial license** — check expiry before long
  calibration runs.
- **External libraries in place.** `external/SparseEcon` (clone, commit
  `9d04a89`, Feb 2023); `external/Moll-codes/` with the scripts downloaded
  from benjaminmoll.com/codes: `HJB_stateconstraint_{explicit,implicit}.m`,
  `huggett_partialeq.m`, `HJB_diffusion_implicit.m`,
  `huggett_diffusion_partialeq.m`, `aiyagari_diffusion_{equilibrium,transition}.m`,
  `lifecycle.m`, plus `HACT_Numerical_Appendix.pdf`.
- **Smoke test passed.** SparseEcon `use_cases/01_huggett/01_baseline/main.m`
  ran end-to-end (~60 s, converged: r = 0.0199, market-clearing residuals
  ~1e-10). Evidence figures in `external/smoke_test_outputs/`.
- **Scaffold done.** `startup.m`, `src/params_default.m` (Step 1 values from
  appendix §8), `src/grids_build.m` (grids only, no operators), stubs for the
  five Step 1 blocks (`income/returns/mortality/bequest/tax_step1.m`, header +
  signature + `error(...)`). Verified: `grids_build(params_default())` gives
  N = 172,200 = 300 × 7 × 82, matching §3.1.

## State after session 2 (2026-06-11, afternoon)

- **Specs moved into the repo**: `docs/implementation_appendix.md`,
  `docs/preliminary_paper_outline.md`.
- **All five Step 1 block modules implemented** (income/returns/mortality/
  bequest/tax `_step1.m`) plus a units layer in `params_default.m`
  (`eur_per_unit`, see discrepancy 8) and pension floor (discrepancy 9).
- **`operator_build.m`**: dispatcher on toggles; static Kronecker blocks
  `Ay` (OU upwind diffusion, reflecting), `Ah` (unit age drift, forward
  upwind; aging out of h_max = certain death with bequest), mortality diag,
  bequest RHS, stationary `p_y` for heirs' entry income. Decomposed so the
  future `transition_solve.m` reuses everything; only the policy-dependent
  wealth-drift block is rebuilt per HJB iteration (in `hjb_solve.m`).
- **`hjb_solve.m`**: implicit upwind scheme, nonuniform a grid,
  state-constraint boundaries, death-with-bequest inhomogeneous term.
  Converges in ~10 iterations on the full 3-D Step 1 model.
- **`kfe_solve.m`**: stationary KFE in point masses m = g .* wx (the upwind
  generator is an exact rate matrix, so conservation is exact); inheritance
  kernel R recycles mortality + h_max outflows into entry at h0 with
  post-tax wealth `max(b - T(b), 0) + G` (linear interpolation onto the a
  grid) and y′ ~ p_y. Solve: fix-one-row on (A′ + R), normalise.
- **`equilibrium.m`**: wrapper grids→ops→HJB→KFE + revenue flow. No tau
  bisection yet.
- **Tests passing** (run via `matlab -batch "startup; test_..."`):
  - `tests/test_aiyagari_limit.m` — age collapsed (Nh = 1), no mortality:
    matches an inline log-OU adaptation of `huggett_diffusion_partialeq.m`
    on identical grids to machine precision (V err ~6e-16, mass err ~4e-15).
  - `tests/test_mass_conservation.m` — full 3-D Step 1 solve (Na = 100):
    exact generator conservation, residual, m ≥ 0, entry inflow = death
    outflow, hump-shaped wealth-by-age profile.

## State after session 3 (2026-06-12)

- **Status-quo equilibrium at full resolution** (Na = 300, N = 172,200):
  ~110 s cold / ~60 s warm, HJB 22 iterations, KFE residual ~6e-17. No
  reordering needed. Figures + `status_quo.mat` in `results/step1/`.
  Levels (theta_b = 1 placeholder): mean wealth 406 k EUR, Gini 0.55,
  top-1% 6.9% — top tail far too thin vs. data, as expected for a
  one-asset single-rate model (see TODO "Model fit").
- **`moments.m`** (top shares via threshold-atom splitting, Gini from the
  Lorenz curve, bwr, revenue, profiles) and **`run_status_quo.m`**.
- **`calibrate_step1.m`**: Illinois root-finder on log10(theta_b) with HJB
  warm starts (`hjb_solve`/`equilibrium` take optional V0) — 3 evals
  instead of ~15 bisections. Converges, BUT the BWR calibration is
  degenerate under constant mortality (see TODO) — theta_b stays at the
  1.0 placeholder; pull Step 2.C forward before calibrating for real.
- **`params_derive.m`**: recomputes model-unit fields from `*_eur` inputs;
  call it after overriding any EUR field in an experiment script.
- **`solve_revenue_balance.m`** (reusable core): INCREMENTAL revenue balance
  (discrepancy 13) Rev(tau*) = Rev_SQ + G·N_entry. **Solves the HJB ONCE and
  sweeps tau with KFE-only re-solves**, because in Step 1 V is invariant to
  tau and G (see next bullet) — fzero (Brent) on tau over cheap KFE solves.
  Cuts the outer loop from ~6 full HJB+KFE solves (~27 min, was the real
  bottleneck — NOT the root-finder) to 1 HJB + ~8 KFE (~2 min). Guards
  against misuse when bequest≠'step1' (net-of-tax motive would break the
  invariance). `run_exp1_grunderbe.m` (G=20k) and
  `run_exp1b_grunderbe_200k.m` (G=200k) are thin wrappers;
  `save_experiment_figures.m` is the shared plotting helper.
  **Exp1 result** (theta_b = 1 placeholder): tau* = 0.3722 (Na=300; 0.3661
  at Na=100). Clean redistribution: mean wealth 406k→401k EUR (−1.4% = the
  MPC channel), bottom-50% 12.7→13.9%, Gini 0.548→0.532, top-1% ~unchanged.
  **Exp1b (G=200k) is INFEASIBLE under revenue balance**: even tau=1 (100%
  above the 400k exemption) raises an increment of 0.061 vs a 0.111 bill.
  The largest revenue-balanceable Grunderbe is ~110k EUR — a headline bound
  (upper bound: ignores donor response, placeholder theta_b + constant
  mortality). The estate-tax base above 400k is too thin to fund a grant
  paid to every entrant beyond ~110k.
- **Step 1 estate tax has zero donor response by construction**: the warm
  glow W(a) values the GROSS estate, and the grant enters only the entry
  kernel, so neither tau nor G enters the HJB — V/policy/generator A are
  policy-invariant. This is what licenses the one-HJB speedup above, and
  it means the whole redistribution runs through the KFE kernel. Note for
  interpretation and for the Step 2 bequest-motive decision (see TODO).
- **MATLAB stability**: a batch run died at a fatal interpreter crash
  (`~/matlab_crash_dump.27552-1`) while the machine was 14.6 GB into swap;
  close the desktop MATLAB during long runs (see TODO).
- **`tests/run_tests.m`** runs the whole suite in one MATLAB session.
- Workflow notes: never pipe long MATLAB batch runs through `tail` (output
  is lost if the process dies); announce runtime estimates for anything
  >10 min (one stationary solve: ~40-60 s at Na = 100, ~2 min at Na = 300).
  Harness background tasks get killed at session handovers (two MATLAB
  runs died this way mid-run): launch anything >5 min detached instead —
  `nohup matlab -batch "..." > results/<log>.txt 2>&1 & disown` — and
  watch the log file. Design long loops so they can resume from partial
  output (e.g. `run_exp1_grunderbe(params, tau_lo, tau_hi)`).

## Conventions adopted

- `external/` is git-ignored (Moll's scripts are loose files, not a repo, so
  submodules were not an option for both); provenance pinned here and in
  README.
- File names lowercase_with_underscores per appendix §6; project root plays
  the role of the appendix's `preliminary_paper/`.
- `params` struct from `params_default.m` is the single source of truth;
  module toggles are strings: `params.income = 'step1'` etc.
- Block-module signatures (consumed by the future `operator_build.m`):
  - `[drift_y, var_y] = income_step1(grids, params)` — OU coefficients on the
    y grid (state kept in **logs**, so no Itô correction needed).
  - `r = returns_step1(grids, params)` — N×1 on the flattened grid.
  - `lambda = mortality_step1(grids, params)` — Nh×1 by age.
  - `W = bequest_step1(grids, params)` — Na×1 warm-glow value.
  - `T = tax_step1(b, params)` — tax due on bequest b (vectorised on values,
    not the grid; used by the inheritance kernel and revenue integration).
- Grids: wealth grid is power-spaced, `a = a_max * u.^a_curv` with
  `a_curv = 2` (refines near a = 0); y grid uniform over μ_y ± 3 stationary
  sd; flattening is column-major over (a, y, h), a fastest — Moll convention.
- Headless runs: `/Applications/MATLAB_R2026a.app/bin/matlab -batch "..."`
  from the project root (runs `startup.m` automatically when started there).
- SparseEcon use cases write run outputs into their own *tracked* `output/`
  directories, so running one dirties the clone. After a smoke run, restore
  with `git checkout -- .` inside `external/SparseEcon` (evidence copies live
  in `external/smoke_test_outputs/`). Never edit anything under `external/`.

## Closest reference script — shortlist and pick

1. **Moll `huggett_diffusion_partialeq.m`** — stationary HJB+KFE joint solve,
   diffusion income, partial equilibrium, borrowing constraint. Exactly our
   Step 1 income/equilibrium blocks, and the natural regression target for
   the "Aiyagari limit" test; lacks the age dimension.
2. **Moll `lifecycle.m`** — life-cycle with log-OU diffusion income, age
   handled by backward time-marching (age = time), terminal value ≈ 0.
   Closest to our dynamics; but no mortality hazard, no bequest motive, no
   stationary KFE with entry, and age-as-time doesn't give the stationary
   age-as-state operator the appendix prescribes.
3. **SparseEcon `use_cases/08_life_cycle/01_one_asset_life_cycle`** — age as
   an actual grid dimension `(a, t)` with OLG birth/death bookkeeping in
   `KF.m`; structurally the right shape. But: adaptive sparse grids (not the
   dense tensor grid of §3.1), two discrete income types (not diffusion),
   and general equilibrium in K.

**Pick: `huggett_diffusion_partialeq.m` as the base skeleton**, extended with
the age dimension as a deterministic unit drift following `lifecycle.m`'s
upwind logic, and SparseEcon's life-cycle `KF.m` as the reference for the
death→entry boundary operator. One line each: (1) is the only candidate with
the full stationary HJB+KFE pattern *and* diffusion income *and* PE, i.e. the
most machinery reusable verbatim; (2) contributes the age treatment but
solves a finite-horizon problem, not a stationary distribution; (3) is the
right economics but the wrong numerical stack — its sparse-grid library
diverges from the appendix's dense Kronecker prescription.

## Next session (development order §7, item 5–7)

1. Status-quo equilibrium at full resolution (Na = 300): run
   `equilibrium(params_default())`, plot wealth density and wealth-by-age
   profile, first sanity-check of levels (§7 item 5). Check runtime; if the
   full-grid sparse solve is slow, consider reordering states (a slowest)
   to cut bandwidth before optimising anything else.
2. `moments.m` (top shares, Gini, bequest-to-wealth ratio) — needed for
   calibration targets.
3. `calibrate_step1.m`: bisection on theta_b to match the target
   bequest-to-wealth ratio (§7 item 6).
4. Then the revenue-balance outer loop and `run_exp1_grunderbe.m` (§7
   item 7).
5. Remaining validation: `test_lifecycle_limit.m` (lambda→0, theta_b→0
   against a deterministic benchmark), no-income-heterogeneity limit.

## Design note: steady states now, transitions later

The dense-Kronecker (Moll) stack was chosen over SparseEcon's adaptive
sparse grids and this does NOT change with transition dynamics in view:
in PE with no time-varying prices, a transition is one backward HJB sweep +
one forward KFE sweep over the same operator (appendix §3.5), i.e. ~2 T
sparse solves of the system we already factor once per HJB iteration —
entirely feasible at N ≈ 1.7e5. Sparse adaptive grids pay off in higher
dimensions or GE fixed points, neither of which we have. What transitions
DO require is already in place: `operator_build` keeps the static blocks
(Ay, Ah, mortality) separate from the policy-dependent wealth-drift block,
and `kfe_solve` exposes the inheritance kernel R, so the forward sweep can
apply (A_t′ + R) at every time step without reassembly.

## Discrepancies vs. appendix

1. **`huggett_diffusion_partialeq.m` runs the OU in *levels*, not logs.** The
   log-OU branch exists but is commented out; the appendix (§2 Block A.1)
   assumes the script's style matches log-OU. Not blocking: we keep the state
   y = log income and exponentiate in the budget constraint, which is cleaner
   than either reference (constant OU coefficients, no Itô term).
2. **N_y = 7 is very coarse for a diffusion.** Moll uses J = 40
   (huggett_diffusion) / J = 15 (lifecycle) points to resolve the
   second-order income operator. With 7 points, dy ≈ 0.63 — the diffusion
   term will be poorly resolved. Kept 7 per §8 for now; cheap to raise, and
   we should revisit before validation runs. (N_y = 7 reads like a leftover
   from a discrete-Markov-chain formulation.)
3. **"Aiyagari limit" regression target (§5.1.3).** Appendix says to compare
   the collapsed model to "Moll's stationary Aiyagari diffusion script", but
   that script (`aiyagari_diffusion_equilibrium.m`) is *general* equilibrium
   (solves for r). Our model is PE with fixed r0, so the correct regression
   target is `huggett_diffusion_partialeq.m` with r set to r0.
4. **"Geometric spacing" of the a grid (§3.1) is ill-defined from a = 0**
   (constant-ratio spacing can't start at zero). Implemented the standard
   HACT alternative: power-law spacing `a_max * u.^curv`, curvature 2, which
   delivers the intended refinement near the constraint. Flagging the
   interpretation, can switch to a shifted-log grid if preferred.
5. **SparseEcon's stack is adaptive sparse grids**, not dense tensor grids +
   Kronecker operators as §3.1 prescribes. So "inherit working solvers from
   SparseEcon" mostly means: reuse its templates/patterns (e.g. `KF.m`
   birth-death bookkeeping, transition-dynamics layout in
   `use_cases/05_transition_dynamics`), while the dense-grid machinery comes
   from the Moll scripts. The appendix's §3 numbers (N ≈ 1.7e5) are all
   dense-grid, consistent with this reading.
6. **Appendix §3.3 transition reference "`use_cases/05_...`"** resolves to
   `use_cases/05_transition_dynamics` — exists, fine.
7. **MATLAB is on a trial license** (R2026a, license "DEMO") — fine for now,
   but worth resolving before the ~100-solve Step 2 calibration runs.
8. **§8 units are internally inconsistent**: `w = 1.0` (normalisation) but
   `a_max`, `F`, `G` in EUR. With w = 1 nobody ever reaches a 4e5 exemption.
   Resolved via explicit `params.eur_per_unit = 4.5e4` (≈ German mean annual
   gross earnings, Step 1 placeholder): EUR-denominated parameters are
   converted to model units in `params_default.m` (`*_eur` fields are the
   inputs, model-unit fields derived). Keeps w = 1 per §8 and V, c well
   scaled. Calibrate `eur_per_unit` properly in Step 2.
9. **Appendix budget has zero income in retirement** (`w e^y 1{h<hR}` only),
   which gives c = 0 and u = -Inf for retirees at the constraint — the HJB
   is ill-posed without a floor. Added `params.pension_eur = 1.35e4` (~30%
   of mean earnings, Grundsicherung-like, Step 1 placeholder). Step 2 should
   replace this with a proper German pension (point system / OECD data).
10. **Heirs' income state at entry is unspecified** (§1.2 only fixes initial
    wealth). Implemented: y′ drawn from the stationary distribution of the
    discretised income process, independent of the parent. Alternative
    (intergenerational y-persistence) would need a copula parameter — Step 2
    decision.
11. **Warm-glow W(a) = θ_b a^(1-γ)/(1-γ) is -Inf at a = 0** for γ ≥ 1.
    Added `params.bequest_shift_eur = 2e4` inside the power (regularisation,
    same role as De Nardi's θ₂). Affects θ_b calibration only marginally.
12. **Terminal age**: aging out at h_max is implemented as certain death
    with bequest (the aging outflow at the last age node feeds the
    inheritance kernel, value side gets W(a)). The appendix says only
    "h_max acts as a hard upper bound"; this is the consistent reading
    (mass must go somewhere), but flagging the choice.
13. **Revenue balance is INCREMENTAL, not total** (user decision
    2026-06-12, overrides appendix §3.4): Rev(tau*) = Rev_status_quo +
    G·N_entry. Rationale: status-quo revenue is a leak out of the
    household sector (no government in the model); §3.4's total balance
    would turn exp1 into a tax cut plus full recycling, mixing
    redistribution with leak elimination (measured: tau* = 0.139 < 0.20,
    mean wealth +5%). Holding the leak constant isolates redistribution.
    Related insight for the paper: with wealth-independent returns (Step
    1), retiming transfers moves aggregate wealth only through MPC
    heterogeneity (compounding cancels — r·W invariant to who holds W);
    with Step 2 r(a), a first-order non-behavioural channel reappears.
