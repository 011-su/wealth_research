# Outline: *Estate-Tax-Financed Capital Grants and the Long-Run German Wealth Distribution*

**Target venue.** DIW Discussion Paper first; journal submission to *Journal of Public Economics* or similar.
**Implementation language.** MATLAB, built on Moll's reference codes (<https://benjaminmoll.com/codes/>) and the SparseEcon toolbox (<https://github.com/schaab-lab/SparseEcon>).
**Companion document.** `implementation_appendix.md` (math + numerics + code structure).

## 1. Scope

A continuous-time heterogeneous-agent life-cycle model with bequests, calibrated to German data, used to compute the long-run distributional effects of four estate-tax / capital-grant policies:

1. ***Grunderbe***: 20,000 EUR at age 18, funded by progressive *Erbschaftsteuer* reform (revenue-balanced in the new steady state).
2. Same revenue-equivalent estate-tax reform, no grant.
3. Revenue-equivalent universal income transfer, paid as a flat lifetime annuity rather than at age 18.
4. Closure of the *Verschonungsregeln* business-wealth exemption in isolation (no other change to the schedule).

For each policy: long-run change in top-1% share, top-10% share, bottom-50% share, Gini, and the wealth-by-age profile. Transition dynamics from pre-reform stationary state. The implementation uses the *two-step strategy* of the appendix: Step 1 (literature-parameterised baseline, adapted directly from Moll reference scripts) for code validation; Step 2 (German-data calibration) for paper results.

## 2. Paper structure

### 1. Introduction (3 pp)

- Erbschaftswelle (cite Tiefensee–Grabka 2017; Westermeier et al. 2016).
- Frühstartrente and the legislative timeline (coalition agreement 2025; BMF *Eckpunktepapier* 2025).
- Research question: long-run effects of estate-tax-financed capital grants on the German wealth distribution.
- Approach: structural continuous-time HA life-cycle model calibrated to German data; four policy scenarios.
- Contribution: first structural quantitative analysis calibrated to German institutional detail (*Steuerklassen*, *Verschonungsregeln*).

### 2. Related literature (1 pp)

HA models with bequests (De Nardi 2004; Benhabib–Bisin 2018); inheritance taxation (Piketty–Saez 2013; Atkinson 2015); *Grunderbe* simulation (Bach 2021; Bach–Sinclair 2026); continuous-time HA framework (Achdou et al. 2022).

### 3. Model (4–5 pp)

**Environment.** State $x = (a, y, h)$ with $a \geq 0$ wealth, $y$ log labor productivity, $h \in [h_0, h_{\max}]$ age. Continuous time. Partial equilibrium ($r$ exogenous).

**Household problem.** Bellman:

$$\rho V(a,y,h) = \max_c \Big\{ u(c) + V_a \cdot \dot a + \mathcal{L}_y V + V_h \Big\} - \lambda(h)\big[V(a,y,h) - W(a)\big]$$

with $u(c) = c^{1-\gamma}/(1-\gamma)$, $\dot a = r(a)\,a + w \exp(y)\,\mathbf{1}\{h < h_R\} + G\,\delta(h - h_G) - c$, generator $\mathcal{L}_y$ from the income process, mortality hazard $\lambda(h)$, warm-glow bequest value $W(a)$.

**Boundary operator.** Inflow at $h = h_0$ equals outflow from mortality. Each new entrant's initial wealth is $a_{\text{init}} = \max(b - T_e(b),\, 0) + G \cdot \mathbf{1}\{\text{grant active}\}$, where $b$ is the bequest received and $T_e(\cdot)$ the estate-tax schedule. Population mass preserved.

**Stationary equilibrium.** Joint density $g(a, y, h)$ solves $\mathcal{A}^\top g = 0$ subject to mass normalisation, where $\mathcal{A}$ is the HJB generator. Aggregate consistency: tax revenue funds grants,

$$\int T_e(b)\, \lambda(h)\, g(a, y, h)\, dx = G \cdot N_{\text{entry}},$$

with $N_{\text{entry}}$ the flow of entrants per unit time.

**Specification choices.** Six parameter blocks; see the implementation appendix for mathematical detail and the two-step structure: preferences, income process, returns, mortality, bequest motive, estate-tax schedule.

### 4. Calibration (3–4 pp)

**Strategy.** Direct calibration where the parameter is observable in external data; indirect inference where it is not (bequest motive in particular). All policy results in the paper use the Step 2 calibration; Step 1 used for sensitivity and validation.

**Targets (Step 2).**

- Aggregate bequest-to-wealth ratio from *Erbschaftsteuerstatistik* adjusted with the Tiefensee–Grabka (2017) upscaling.
- Wealth-by-age profile from SOEP and Bundesbank PHF.
- Top-decile wealth share from Albers–Bartels–Schularick (2024).
- Bequest-share by recipient income decile from Westermeier et al. (2016).

**Calibration table.** [Table 1: each parameter, Step 1 value, Step 2 value, source, target moment if indirect.]

### 5. Status quo (2–3 pp)

Validation: stationary $g(a, y, h)$ under current ErbStG. Compare to:

- ABS top-share series (out-of-calibration: top-1%, top-5%);
- SOEP / PHF wealth-by-age profile (also a direct target, so this is a consistency check);
- Bundesbank aggregate household wealth (level check after scaling).

[Figure 1: status-quo wealth density. Figure 2: wealth-by-age profile, model vs. data.]

### 6. Policy experiments (6–7 pp)

For each of the four scenarios:

1. State the policy change (schedule, parameters).
2. Solve for the new stationary equilibrium under revenue balance (where applicable).
3. Compute transition dynamics from pre-reform $g$ to new stationary $g$.
4. Report: top shares, bottom-50% share, Gini, wealth-by-age profile, transition half-life.

**6.1 Grunderbe (Experiment 1).** $G = 20{,}000$ EUR at $h = 18$. Revenue from estate-tax reform: progressive schedule replacing the current *Steuerklasse*-specific scheme, with rate and exemption structure calibrated so that $\int T_e \lambda g = G \cdot N_{\text{entry}}$ in the new steady state. Outer-loop search for the revenue-balanced schedule.

**6.2 Revenue-equivalent estate-tax reform (Experiment 2).** Same estate-tax reform as 6.1, but $G = 0$. Revenue rebated lump-sum to all decedents' heirs proportional to bequests (mechanical, no grant). Isolates the redistribution-through-taxation channel.

**6.3 Revenue-equivalent universal income transfer (Experiment 3).** Same total revenue paid as a constant flow $\tilde G$ to all working-age households. Isolates the timing-of-transfer channel (mid-life vs. age 18).

**6.4 Verschonungsregeln closure (Experiment 4).** ErbStG unchanged except the business-wealth reduction is set to zero. No grant. Reports the marginal contribution of the *Verschonungsregeln* margin alone.

**Summary results.** [Table 2: top-1%, top-10%, bottom-50%, Gini, transition half-life, across all four experiments. Figure 3: stationary densities overlaid. Figure 4: transition paths of top-1% share.]

### 7. Discussion and robustness (2 pp)

- Sensitivity to the bequest specification: warm-glow vs. De Nardi luxury form.
- Sensitivity to the returns specification: heterogeneous (Step 2 default) vs. single safe rate (Step 1).
- Behavioural margins absent from the model: savings response of donors, business-form conversion around *Verschonungsregeln* thresholds. Discuss qualitatively.
- Comparison to Bach (2021) reduced-form simulation: where the structural model agrees and where it disagrees.

### 8. Conclusion (1 pp)

## Appendices to the paper

- **A. Data.** Source-by-source description of German calibration data (Erbschaftsteuerstatistik, SOEP, PHF, ABS, Destatis): definitions, harmonisation, vintage.
- **B. Computational details.** Distilled from `implementation_appendix.md`.
- **C. Additional results.** Sensitivity figures and tables.

## 3. Deliverables

1. MATLAB code with Step 1 implementation, adapted from Moll/SparseEcon reference scripts, passing all validation checks (see appendix Validation section).
2. Step 1 results: status-quo equilibrium and four policy experiments. Quick sanity check.
3. MATLAB code with Step 2 implementation passing extended validation.
4. Step 2 results: paper Tables 1–2 and Figures 1–4, plus sensitivity tables.
5. Reproducibility material: parameter files, scripts, README. Released as a GitHub repository attached to the DIW Discussion Paper.
