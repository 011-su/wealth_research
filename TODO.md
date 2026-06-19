# TODO — outstanding tasks

Working list of known gaps and decisions deferred. Numbers in brackets refer
to the "Discrepancies vs. appendix" section in `CLAUDE.md`.

## Spec overhaul (2026-06-18) — the appendix now describes a different model

The user replaced `docs/implementation_appendix.md` with a much richer spec
(4-D state `(a, z_p, z_ε, h)`, life-cycle ψ(h), full ErbStG with δ-shift,
fertility-based inheritance, indirect inference, ages 0–78, no two-step
staging). The current code implements a Step-1-like subset. Bringing the code
up to the new contract is now the main roadmap; in rough dependency order:

- [ ] Income: split into persistent (Rouwenhorst, N_p=7) + transitory
  (Tauchen, N_ε=5) components, generators via `logm` of annual transition
  matrices, Kronecker-sum Λ_y (§A.1.1–3). Watch the transitory generator
  validity (near rank-1; use the (P−I)/Δt fallback, §A.1.2).
- [ ] Life-cycle earnings profile ψ(h) quartic from FSS (§A.1.5); age grid
  0–78, retirement collapse at h_ret=43 to replacement rate 0.55 (§A.1.7).
- [ ] Bequest motive: De Nardi W(a)=φ(a−â)^(1−γ)/(1−γ), calibrate (φ,â) by
  indirect inference (`fminsearch`, §F) on BWR / wealth-share-65+ / below-
  exemption-share. (Supersedes the degenerate one-param BWR calibration.)
- [ ] Full ErbStG Steuerklasse-I schedule with δ rate-shift, Verschonungs-
  regeln β(a), €400k→€200k exemption (§D); outer-loop bisection on δ.
- [ ] Inheritance kernel with heir-count f_n(n|h), childless pool, Young
  lottery (§C); fertility data from Destatis.
- [ ] Re-map experiments to the new numbering (exp2 = retirement transfer,
  exp3 = annual transfer, exp4 = Verschonung closure with endogenous G).
- [ ] Decide whether to keep the `params.*='step1'|'step2'` toggle scaffold
  or refactor to the single-target structure the new spec implies.

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
- [ ] **Warm glow over gross vs. net bequest.** Step 1's W(a) values the
  gross estate, so the estate tax and the grant do not enter the HJB at
  all: zero donor savings response by construction (observed: V is
  policy-invariant, warm-started solves converge in 1 iteration). If the
  paper wants a donor margin, Step 2's De Nardi form should value the
  net-of-tax bequest W(a − T_e(a)) — decide alongside Block D.2 and flag
  in the discussion section either way (outline §7 lists the donor
  response as a qualitative margin; in Step 1 it is exactly zero).

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

## Memory / environment

- [ ] **Full-res (Na=300) direct solves exhaust the 16 GB Mac** — root
  cause of the repeated MATLAB crashes (5+ since Jun 12, all "Trace trap"
  under heavy swap), NOT a code bug. The sparse LU fill-in of the (a,y,h)
  Kronecker system at N≈1.7e5 needs more RAM than is free. Mitigations:
  (a) run experiments at Na≤200; (b) reboot to clear the swap backlog
  before a full-res run; (c) reorder states (a slowest) to cut fill-in
  bandwidth; (d) longer term, more RAM or an iterative KFE solve. Resolved
  in practice by moving full-res runs to the HU pool machine (62 GB).
- [ ] **KFE direct solve is slow at full res (~10 min/solve at Na=300)** —
  the inheritance kernel R couples every dying state to the h0 entry block,
  destroying the operator's band structure, so `(A'+R)\rhs` has heavy LU
  fill-in. Fine for one solve, but the revenue-balance sweep (fzero, ~10
  KFE) then takes ~2 h, and the new-spec model (N≈1.4e6) would be
  intractable this way. **Fix: iterative KFE solve** — `bicgstab`/`gmres`
  on (A'+R)' with an ILU or block preconditioner, or reorder states (a
  slowest) so R's coupling stays narrow. Appendix §A.1.9 (annotated) and
  the build plan already flag bicgstab. Do before scaling to the full model.
- [ ] **HU pool machine kills jobs silent on stdout for ~28 min** (idle
  reaper; not OOM, not ulimit — see memory note `remote-server-access`).
  Worked around for the revenue-balance sweep by printing after every KFE
  solve (commit 3cf356b). Any long silent computation launched there must
  emit a periodic heartbeat, or run under a wrapper that does.

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

- [x] **DECIDED 2026-06-12: incremental revenue balance** Rev(tau*) =
  Rev_status_quo + G·N_entry (leak held constant across scenarios), not
  the appendix §3.4 total balance — see CLAUDE.md discrepancy 13.
  Implemented in `run_exp1_grunderbe.m`; apply the same convention to
  exp2–exp4 when they are written.
- [ ] Timing-channel note for the paper (from the same discussion): with
  wealth-independent returns, retiming transfers affects aggregate wealth
  only via MPC heterogeneity (compounding cancels in the aggregate); with
  Step 2 r(a) heterogeneity a non-behavioural channel reappears. Make
  this explicit when interpreting exp1 vs exp3.

- [x] Revenue-balance outer loop: `solve_revenue_balance.m` (incremental
  condition, one-HJB speedup, fzero on tau). Wrappers: exp1 (20k), exp1b
  (200k). exp2–exp4 still to write (reuse the core).
- [ ] **exp2–exp4 must re-solve the HJB per parameter** once their policies
  enter the HJB (e.g. universal transfer changes the budget; a net-of-tax
  bequest motive). The one-HJB shortcut in `solve_revenue_balance.m` is
  specific to Step 1's gross-estate warm glow + entry-kernel grant; it
  guards with a warning when bequest≠'step1' but the assumption should be
  rechecked per experiment.
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
