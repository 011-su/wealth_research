# Qualitative results & key insights

Living record of the main findings from the WP0 model runs. Numbers are
**preliminary**: the model is the Step-1 machinery with Step 2.C (Destatis
age-rising mortality) pulled forward, and the bequest-motive strength
`theta_b` is still at the placeholder value 1.0 (not yet calibrated to data).
Magnitudes will move with calibration; the *directions and mechanisms* below
are robust unless flagged.

## Setup

- Continuous-time HA life-cycle model, partial equilibrium, German units
  (1 model unit = €45k mean earnings).
- **Experiment 1 (Grunderbe):** a lump-sum capital grant `G` to every entrant
  at age 18, funded by a **flat surtax `tau_add` on *all* bequests (no
  exemption)**, on top of the unchanged status-quo estate tax (held constant
  as a leak). Revenue balance: `tau_add · E[bequest·death rate] = G · N_entry`.
- "Step 2" below = **Destatis 2022/24 age-rising mortality** (mean age at
  death ≈ 81), replacing the constant-hazard Step-1 baseline (which put ~20%
  of "deaths" at the forced terminal age 99).

## Key quantitative results — Step 2 mortality

Status quo (G = 0): bottom-50% wealth share **13.8%**, Gini **0.532**,
mean wealth **€478k**, top-10% 36.4%, top-1% 6.5%.

| Grunderbe | required surtax `tau_add` | bottom-50% | Gini | mean wealth |
|---|---|---|---|---|
| €0 (SQ) | – | 13.8% | 0.532 | €478k |
| €20k | 7.3% | 14.3% | 0.525 | €476k |
| €100k | 36.8% | 16.4% | 0.500 | €469k |
| €200k | 74.3% | 18.6% | 0.480 | €465k |

Figure: `results/step2/step2_comparison.png` (wealth-over-age and the
non-cumulative wealth density, overlaid by grant size). Step 1 (constant
mortality) funds the same grant at a *lower* surtax (€100k: 30.9%,
€200k: 62.6%) — see Insight 4.

## Insight 1 — The grants are progressive, via the bottom half

Each grant **raises the bottom-50% share and lowers the Gini**, while the
**top-10% and top-1% shares barely move**. Redistribution operates by lifting
the asset-poor bottom half, *not* by compressing the top. At €200k the
bottom-50% share rises ~5 pp (13.8 → 18.6%) and the Gini falls ~0.05
(0.532 → 0.480). The effect scales steadily with grant size.

## Insight 2 — What the wealth distribution shows

(Right panel of the figure.) The grant **moves households off the borrowing
constraint**: the density spike at `a ≈ 0` shrinks as `G` rises, and mass
piles into a **hump at the grant level** (clear at €100k, pronounced at
€200k where the low-wealth spike all but disappears — every entrant now starts
with ~€200k). The upper tail is essentially unchanged. So the grant compresses
the **cross-sectional** distribution from below.

(Left panel.) The **wealth-over-age profiles barely move** across grant sizes —
the hump-shaped accumulate-to-retirement-then-decumulate life-cycle shape is
preserved. The grant changes who holds wealth in the cross-section, not the
average life-cycle trajectory.

## Insight 3 — The €400k exemption, not the inheritance flow, caps the grant

The model's bequest flow is realistic (mean bequest ≈ €330k per decedent,
aggregate ≈ €340 bn/yr, ~1.0M deaths/yr — in line with German data). An estate
tax *above the €400k exemption* nonetheless funds only a small grant (ceiling
≈ €84k even at a 100% top rate) because the exemption sits near the **72nd
percentile of estates** and the model's **top tail is too thin** (top-1% 6.5%
vs ~27% in data) — so the taxable base above €400k is narrow. Removing the
exemption (taxing *all* bequests) widens the base ~2.5× and makes even a €200k
grant fundable.

## Insight 4 — Realistic mortality raises the funding cost

Under Step 2 mortality the **required surtax is higher** for every grant
(€100k: 36.8% vs 30.9%; €200k: 74.3% vs 62.6%). Mechanism: households die
**old** (~81) after **decumulating through a long retirement**, so they die
*poorer than their peak* and the annual bequest flow is smaller
(bequest-to-wealth ratio 0.0092 vs 0.0205 under constant mortality), shrinking
the surtax base. (This same decumulation is why the inheritance flow looks
"small" under realistic mortality, and it is what makes the `theta_b`
bequest-motive calibration well-identified.)

## Insight 5 — Feasibility frontier and political magnitude

A €200k universal Grunderbe **is** fundable from bequests alone, but only with
an **extreme** flat rate on *every* inheritance: ~63% (Step 1) to ~74%
(Step 2). Under realistic mortality the rate is near the hard cap (surtax plus
the status-quo schedule reaches a 100% top marginal estate rate at
`tau_add` = 0.8), so **€200k is roughly the largest grant the broad-based
bequest surtax can finance.** Smaller grants are comfortably feasible (€20k:
7%, €100k: 37%).

## Known artifacts / open issues (to address later)

- **Aggregate (mean) wealth declines with the grant** (−0.4% at €20k to −2.7%
  at €200k, Step 2). This is an important model feature but is **treated as an
  artifact for now** — its sign and size are not yet trustworthy because they
  depend on the **uncalibrated bequest motive** (`theta_b`) and on the
  partial-equilibrium, wealth-independent-returns structure. *Candidate*
  mechanism to verify after calibration: the grant moves resources to
  high-MPC young households who consume more, lowering the steady-state stock
  (with constant returns the compounding channel cancels, so any effect is the
  MPC/consumption response). Revisit once `theta_b` is calibrated and Step 2
  returns/GE are in. See `TODO.md`.
- Top tail far too thin vs data (top-1% 6.5% vs ~27%); fixed in plan by Step 2
  heterogeneous returns + De Nardi bequests.

## Caveats / not yet calibrated

- `theta_b` is a placeholder; calibrating it (now well-identified under
  Step 2.C mortality) will shift the bequest flow, the surtax base, and the
  rates above.
- Step 2 income (persistent+transitory, life-cycle profile), heterogeneous
  returns, full ErbStG schedule, De Nardi bequest form: not yet implemented
  (see `TODO.md`, "Spec overhaul").
- No behavioural donor response in Step 1 (warm glow values the gross estate);
  partial equilibrium (fixed `r`).
