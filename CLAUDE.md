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

## Appendix overhaul (2026-06-18) and code divergence

`docs/implementation_appendix.md` was **replaced wholesale** by the user with a
much more detailed spec ("Technical Appendix for Implementation", §A–§J, EUR
references retitled). It is no longer the model the code implements:

- **4-D state** `(a, z_p, z_ε, h)`: income split into a persistent OU
  (Rouwenhorst, N_p=7) and a transitory OU (Tauchen, N_ε=5), generators via
  matrix-log of annual transition matrices. Code has a single OU `y` (N_y=7).
- **Ages 0–78** (biological 22+), retirement at h_ret=43 via collapse to a
  single deterministic state at replacement rate 0.55. Code: h0=18, h_max=99,
  hR=65, Nh=82, pension floor.
- **Life-cycle earnings profile** ψ(h) (quartic, FSS). Code: none.
- **Full ErbStG Steuerklasse-I** progressive schedule with a δ rate-shift,
  Verschonungsregeln, €400k→€200k exemption. Code: flat `tau0·max(b−F,0)`.
- **Fertility-based inheritance**: heir-count f_n(n|h), childless pool,
  Young-lottery kernel. Code: single-heir split, no fertility.
- **Indirect inference** for (φ, â) bequest params (De Nardi form). Code:
  one-param warm glow, degenerate BWR calibration.
- **Two-step (Step 1 / Step 2) staging dropped.** The code is built around
  `params.income='step1'|'step2'` toggles; the new spec is single-target.
- **Experiments renumbered/changed**: exp2 is now a retirement transfer (not
  "reform, no grant"); exp1a/1b are the G=20k/200k Grunderbe.

So the current code implements a **Step-1-like subset** of the new spec.
Closing the gap is now the main roadmap item — see TODO "Spec overhaul".
**Three of our run's insights were folded INTO the new appendix**: incremental
revenue balance (§G), the one-HJB-solve optimisation (§B.4 "subtle point on
coupling", §E.2), and the G=200k infeasibility (§G.7, test 12).

## Discrepancies vs. appendix (reconciled against the 2026-06-18 appendix)

Status tags: RESOLVED (new appendix specifies it), MOOT (model changed),
OPEN (still live), ADOPTED (folded into the spec). Section refs are to the
NEW appendix.

1. OU levels-vs-logs — **MOOT.** New spec is log-productivity with an explicit
   Jensen correction (§A.1.6); the reference-script quirk no longer bites.
2. N_y=7 coarse — **SUPERSEDED.** Income is now 35 states (7×5) via
   Rouwenhorst+Tauchen+matrix-log (§A.1.1–3); Rouwenhorst handles persistence
   at low N. New risk instead: matrix-log of the near-rank-1 *transitory*
   matrix often has no valid generator — **appendix patched (§A.1.2)** to apply
   the §A.1.1 validity check + (P−I)/Δt fallback there too.
3. Aiyagari-limit target — **MOOT.** New validation is an Achdou Huggett
   replication (test 1, §J), not the GE Aiyagari diffusion script.
4. Geometric spacing — **RESOLVED.** New §A.1 gives an explicit
   exponential-stretch grid (λ=5), well-defined from a_min=0. Our power-law
   grid is a valid alternative; can switch to the appendix formula.
5. SparseEcon sparse-vs-dense — **PARTLY OPEN.** New spec still prescribes a
   dense tensor grid (1.38M states, §A.1.9) on SparseEcon `lib/` primitives,
   not adaptive grids — same reading, now compounded by a real memory wall:
   **appendix patched (§A.1.9)** with a memory/runtime warning from our crashes.
6. Transition ref `use_cases/05` — **MOOT/fine** (still §E, §I).
7. MATLAB trial license — **OPEN** (environmental); now also the full-model
   memory wall (TODO "Memory / environment").
8. Units EUR-vs-w=1 — **REAL, now patched.** New spec mixes model-unit wealth
   (units of mean earnings) with EUR tax thresholds, conversion left undefined.
   **Appendix patched (§A.1.6)** defining ȳ_EUR (≈€45k); our
   `eur_per_unit=4.5e4` implements exactly this.
9. Zero retirement income — **RESOLVED.** New §A.1.7 gives retirement income =
   0.55 replacement rate; HJB well-posed. Our pension floor is superseded;
   code should adopt the replacement-rate spec.
10. Heirs' income at entry — **RESOLVED.** New §C.1/§C.4: heirs draw z_p, z_ε
    from the ergodic distributions (benchmark), with a persistent-transmission
    sensitivity (IG elasticity 0.32). Matches our stationary-draw choice.
11. Warm-glow −∞ at a=0 — **RESOLVED.** New §A.2 W(a)=φ(a−â)^(1−γ)/(1−γ) for
    a≥â + floor below; â is calibrated (§F). Our `bequest_shift` is the analog,
    but the motive is now the calibrated De Nardi form, not one-param warm glow.
12. Terminal age — **OPEN (minor).** New §A.2 sets V(·,H)=W at H=78 but doesn't
    fully specify the KF mass handling at H; our forced-death-with-bequest
    reading is consistent. Left as-is.
13. Revenue balance incremental — **ADOPTED.** New §G *is* our incremental
    convention (T_e^SQ held constant, only T_e^grant funds the grant). No
    longer a discrepancy — it's the contract. Paper note still stands: with
    wealth-independent returns retiming moves aggregate wealth only via MPC
    heterogeneity (compounding cancels); with heterogeneous r(a) a first-order
    non-behavioural channel reappears.

**Appendix edits made this session** (2026-06-18): §A.1.2 transitory-generator
validity/fallback; §A.1.6 EUR units bridge; §A.1.9 memory/runtime realism;
§B.4 δ_max reconciled to 0.70 (was "0.4", contradicting §G.7's hard cap, with
0.40 kept as a political reference band).
