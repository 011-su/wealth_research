# CLAUDE.md — WP0 preliminary paper implementation

Continuous-time HA life-cycle model with bequests, German estate-tax /
capital-grant experiments. Specification: `implementation_appendix.md` +
`preliminary_paper_outline.md` (the appendix is the contract; both currently
in `~/Downloads/files(8)/`). MATLAB, building on Moll reference codes and
SparseEcon (both in `external/`, read-only, git-ignored).

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

## Next session (development order §7, item 2–3)

1. Run `huggett_diffusion_partialeq.m` logic through the new scaffold
   (params/grids) unchanged — confirm nothing broken (§7 item 2).
2. Start `operator_build.m`: implement `income_step1.m`, `returns_step1.m`,
   build the (a, y) operator on the dense grid, reproduce the Huggett
   stationary result as a regression test.
3. Then add the age dimension (§7 item 3): forward upwind on h, no-mortality
   no-bequest limit, check wealth-by-age profile.

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
