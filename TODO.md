# TODO — outstanding tasks

Working list of known gaps and decisions deferred. Numbers in brackets refer
to the "Discrepancies vs. appendix" section in `CLAUDE.md`.

## Model realism

- [ ] **Pensions.** Replace the flat Grundsicherung-like floor
  (`params.pension_eur = 1.35e4`, Step 1 placeholder added because the
  appendix budget has zero retirement income, which makes the HJB ill-posed
  [9]) with a realistic German pension: earnings-related via the point
  system (pension ∝ average lifetime earnings ≈ function of y at
  retirement), parameterised from OECD *Pensions at a Glance* / DRV data
  (`data/oecd_pensions_at_a_glance_germany.pdf`). Likely needs pension
  claims as a function of the income state at hR, not a flat amount.
  Belongs to Step 2 Block A (income), but is its own decision.
- [ ] **Heirs' income at entry** is currently an independent draw from the
  stationary income distribution [10]. Add intergenerational persistence
  (copula or AR coefficient from SOEP father-son earnings elasticity ~0.3?)
  in Step 2.
- [ ] **Units.** `eur_per_unit = 4.5e4` is a placeholder [8]; calibrate to
  actual German mean gross earnings (Destatis/SOEP vintage to match other
  targets) in Step 2.
- [ ] **Inter-vivos gifts** are absent: the model's bequest flow at death is
  compared to data that partly measures gifts (Erbschaftsteuerstatistik
  covers Schenkungen too). Decide treatment when fixing the Step 2
  calibration targets.

## Calibration

- [ ] **theta_b target.** `target_bwr` must be pinned down from
  Erbschaftsteuerstatistik + Tiefensee–Grabka (2017) upscaling before
  Step 1 results are quotable (~3% per year).
- [ ] **theta_b is weakly identified under Step 1 mortality — and the BWR
  calibration is degenerate.** With the constant hazard, decedents sample
  the population almost uniformly, so bwr ≈ lambda_bar regardless of
  theta_b (measured: 0.0204–0.0255 for theta_b in 1e-3…1e3). Hitting even
  a feasible target (0.024) forces theta_b ≈ 327, which is pathological:
  mean wealth jumps to 1.45 M EUR and the Gini collapses to 0.25 (data:
  ~0.3–0.4 M EUR and ~0.77). Decision 2026-06-12: keep theta_b = 1.0
  placeholder as the default for all Step 1 machinery validation;
  `calibrate_step1.m` works and is tested, but a meaningful calibration
  **requires pulling Step 2.C forward (Destatis Sterbetafel spline, §7
  item 10, the smallest Step 2 change)** so that decedents are old and
  wealth-rich and the BWR becomes informative about the bequest motive.
  Do this BEFORE quoting any experiment numbers.
- [ ] Step 2 calibration targets (`data/targets_step2.mat`): bequest-share
  by recipient decile (Westermeier et al. 2016), wealth share 65+,
  ABS top-share validation series.

## Numerics

- [ ] **N_y = 7 is very coarse for the income diffusion** [2]; raise (e.g.
  15–40) and check stability of results before any quotable run.
- [ ] If the full-grid solve becomes a bottleneck: reorder states (a
  slowest) to cut the kernel/aging bandwidth, or factor the HJB LHS once
  per iteration block. Profile first.
- [ ] Wealth grid: power spacing with curvature 2 chosen [4]; sensitivity
  check (curvature, a_max) once moments are stable.
- [ ] Small mass spike at the a_max boundary (reflecting-top artifact,
  visible in the status-quo density at ~5e6 EUR). Harmless at current
  levels (~1e-5 density) but re-check after Step 2 fattens the top tail;
  raise a_max if it grows.

## Model fit (known Step 1 limitations, for the paper's discussion)

- [ ] Top tail far too thin vs. German data: status quo (placeholder
  theta_b = 1) gives top-1% ≈ 7%, Gini ≈ 0.55 vs. ~27% and ~0.77 in ABS
  data. Expected for a one-asset model with a single safe rate; the Step 2
  blocks (heterogeneous returns, De Nardi bequests) are the planned fix —
  document how far they close the gap.

## Validation (appendix §5.1, remaining)

- [ ] `test_lifecycle_limit.m`: lambda→0, theta_b→0 against a deterministic
  life-cycle benchmark (§5.1.4).
- [ ] No-income-heterogeneity limit sigma_y = 0 (§5.1.5) — needs a guard
  for the degenerate y generator.
- [ ] `test_step1_reproducible_from_step2.m` once Step 2 modules exist.

## Experiments / infrastructure

- [ ] Revenue-balance outer loop (bisection on tau) in `equilibrium.m` or
  experiment scripts (§3.4); then `run_exp1_grunderbe.m` … exp4.
- [ ] `transition_solve.m` (backward HJB + forward KFE over the cached
  operator blocks; design already supports it — see CLAUDE.md design note).
- [ ] **MATLAB trial license** [7]: resolve before Step 2 calibration
  (~100 solves) or long transition runs.
- [ ] **Investigate MATLAB fatal crash of 2026-06-12 17:37**
  (`~/matlab_crash_dump.27552-1`: trace trap in the interpreter thread
  during a forced shutdown, mid-batch-run). Suspects: memory pressure
  (long-running desktop MATLAB + batch solves needing several GB for the
  sparse LU at N = 1.7e5) or trial-license session limits. The same run
  had slowed ~3x before crashing, consistent with swapping. Mitigation
  meanwhile: close the desktop MATLAB during long batch runs; keep runs
  detached and resumable.
