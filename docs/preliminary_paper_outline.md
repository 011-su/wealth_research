# Preliminary Working Paper Outline

## *Estate-Tax-Financed Capital Grants and the Long-Run German Wealth Distribution: A Quantitative Analysis*

**Working title.** *Grunderbe, Estate Taxation, and the German Wealth Distribution: A Quantitative Analysis in a Calibrated Heterogeneous-Agent Model.*

**Target venue.** First as a DIW Discussion Paper. Subsequent journal submission to the *Journal of Public Economics*, *European Economic Review*, *Quantitative Economics*, or *Review of Economic Dynamics*.

**Purpose for the proposal.** Demonstrate three things concretely: (i) the candidate can implement, calibrate, and analyze a continuous-time heterogeneous-agent model in the Achdou–Han–Lasry–Lions–Moll tradition; (ii) the candidate engages substantively with the German policy debate around the *Grunderbe*; (iii) the technical infrastructure for WP2 (HJB–KF solver, transition algorithm, policy counterfactuals) is built and validated before the funded period begins. The paper is self-contained — it uses literature calibration values rather than the WP1 filter — and therefore has independent value.

> **Companion document.** The full discretization, algorithmic specification, inheritance-kernel equation, tax function, module-to-codebase mapping, and unit-test plan are in *Technical Appendix for Implementation* (separate file). Section references of the form "App. §X" below point to that document.

---

## 1. Motivation and contribution

The German wealth distribution is among the most concentrated in the OECD, with the top one percent holding a share comparable to that of the United States and far above the euro-area average. Proposals to broaden wealth ownership through a universal capital grant — the *Grunderbe* of around €20,000 paid at age eighteen, financed by reform of the *Erbschaftsteuer* — have been advanced by Bach (2021) and the surrounding DIW policy work, and have entered the political programmes of the SPD and the Greens in different forms. Yet the existing analyses of these proposals rely on accounting decompositions and partial-equilibrium reasoning. There is no published quantitative analysis of the *Grunderbe* in a structural heterogeneous-agent model calibrated to German data.

This paper provides one. Using a continuous-time life-cycle heterogeneous-agent model with heterogeneous capital returns and warm-glow bequest motives, calibrated to German data, I compute the long-run distributional effects of four redistribution policies. **A central design constraint, motivated by the policy logic of the proposal: the grant is funded entirely by additional *Erbschaftsteuer* revenue raised over and above the status-quo level — no income-tax recourse, no general-fund subsidy. Aggregate estate-tax revenue under the reform decomposes as $T_e = T_e^{\mathrm{SQ}} + T_e^{\mathrm{grant}}$: the first term continues to flow to general government at its status-quo level and is held constant across scenarios for comparability, while only the second funds the grant.** Two grant-size targets are tested — $G = €20{,}000$ (the Bach 2021 headline) and $G = €200{,}000$ (an order-of-magnitude larger stress test) — with the tax-reform parameters adjusted in an outer loop to deliver the corresponding additional revenue $T_e^{\mathrm{grant}} = G \cdot N_{h_0}$ in each case. The four policies are: (i) the *Grunderbe* itself — progressive *Erbschaftsteuer* reform with the additional revenue paid as an age-18 grant; (ii) the same *Erbschaftsteuer* reform with the additional revenue rebated as a lump-sum transfer at retirement instead, isolating the timing channel; (iii) the same reform with the additional revenue rebated as a uniform annual income transfer over working life, isolating the lump-sum-vs-flow channel; (iv) closure of the *Verschonungsregeln* business-wealth exemption alone, with the additional revenue paid as an age-18 grant (the only experiment without a target $G$ — the grant size is whatever Verschonungsregeln closure raises endogenously).

For each policy I report effects on the wealth distribution (top 1%, top 10%, bottom 50% shares, Gini), on cross-sectional inequality decomposed by age, and on intergenerational mobility, with transition dynamics computed from the pre-reform stationary distribution.

The contribution is threefold. First, the paper provides the first structural quantitative analysis of the *Grunderbe* in a model with realistic German features. Second, it isolates the contribution of the *Verschonungsregeln* margin, which the policy debate often treats jointly with rate reform. Third, it characterizes the transition dynamics — how long the reform takes to bind on the steady state — which is crucial for policy design but absent from the existing German discussion.

---

## 2. Model

The model is a continuous-time heterogeneous-agent life-cycle economy in the Achdou et al. (2022) tradition. Households are heterogeneous in financial wealth $a$, a persistent labor-productivity component $z_p$, a transitory labor-productivity component $z_\varepsilon$, and age $h \in [0, H]$. There is no aggregate uncertainty.

**Household problem.** A household of age $h$ with wealth $a$ and productivity state $(z_p, z_\varepsilon)$ chooses consumption $c$ to solve

$$\rho V = \max_c \Big\{ u(c) + \partial_a V \cdot s + \tfrac{1}{2}\sigma_r^2 a^2 \partial_{aa} V + \mathcal{L}_{z_p} V + \mathcal{L}_{z_\varepsilon} V + \partial_h V \Big\} - m(h)\big[V - W(a)\big],$$

with savings rate $s = r(a)\, a + w\, y(z_p, z_\varepsilon, h) - c$, CRRA utility $u(c) = c^{1-\gamma}/(1-\gamma)$, wealth-dependent expected return $r(a)$, idiosyncratic return volatility $\sigma_r$, generators $\mathcal{L}_{z_p}$ and $\mathcal{L}_{z_\varepsilon}$ for the two productivity components, age-specific mortality hazard $m(h)$, and warm-glow continuation value $W(a)$ at death. Boundary conditions (state-constraint at $a = a_{\min}$, zero-flux at $a = a_{\max}$, terminal $V(\cdot, H) = W(a)$) are spelled out in App. §A.2.

**Labor productivity and earnings.** Following Fuchs-Schündeln, Krueger & Sommer (2010, *Review of Economic Dynamics*), log household earnings decompose into a deterministic life-cycle profile $\psi(h)$, a persistent stochastic component $z_p$, and a transitory iid component $z_\varepsilon$:

$$\log y(z_p, z_\varepsilon, h) = \psi(h) + z_p + z_\varepsilon - \tfrac{1}{2}(\sigma_{z_p}^2 + \sigma_{z_\varepsilon}^2)$$

(the Jensen correction normalizes $\mathbb{E}[y] = \exp(\overline{\psi(h)})$ on average over the working life, with $\psi$ itself rescaled so that the population mean of $y$ equals one). $\psi(h)$ is a quartic in age following the FSS specification. Earnings are $w \cdot y$ where $w$ is the wage per efficiency unit of labor — a units constant since the model has no labor-market clearing.

FSS specify the persistent component as a random walk and the transitory as iid. For embedding in a stationary continuous-time HACT model we approximate both as Ornstein–Uhlenbeck processes: the persistent component with very slow mean reversion (effectively a random walk over the working life) and the transitory with very fast mean reversion (effectively iid each year). The OU parameters are pinned by the FSS innovation variances and a small mean-reversion correction; details in App. §A.1 and §C.4.

**Warm-glow bequest motive.** The value of leaving wealth $a$ at death is $W(a) = \phi \cdot (a - \underline{a})^{1-\gamma}/(1-\gamma)$ for $a \geq \underline{a}$, with $\phi$ governing the strength of the bequest motive and $\underline{a}$ a threshold below which the motive is inactive (De Nardi, 2004). $W(a)$ is the value to the *decedent* of leaving $a$; the actual after-tax wealth that lands on heirs is set by the tax-and-division operator $\mathcal{B}$ below.

**Heterogeneous returns.** Following Fagereng et al. (2020), the expected return is wealth-dependent: $r(a) = r_0 + \kappa \cdot \mathbb{1}\{a \geq a^*\} + \chi (a/\bar a)^\eta$, with baseline return $r_0$ (calibrated to 0.02), a premium $\kappa$ above a wealth threshold $a^*$ capturing access to private-equity / business returns, and a smooth scaling term. Idiosyncratic return shocks of variance $\sigma_r^2$ enter both the wealth diffusion in the HJB and the corresponding term in the KF.

**Demographics.** Households enter at age $h = 0$ (biological age 22) with initial wealth $a_0$ determined by the inheritance kernel, age deterministically at unit rate, and die with hazard $m(h)$ calibrated to German life tables (Destatis). Total population is stationary by construction (births balance deaths).

**Retirement.** At model age $h_\mathrm{ret} = 43$ (biological 65), labor productivity becomes deterministic and frozen at $y_\mathrm{ret} = \mathrm{repl} \cdot \mathbb{E}[y]_{\mathrm{working}} = \mathrm{repl}$ (since mean is normalized to one). $\mathrm{repl} = 0.55$ baseline (OECD net replacement rate for Germany). Workers transition to retirement deterministically at age $h_\mathrm{ret}$, with all stochastic productivity components absorbed into the deterministic pension level.

**Government.** The government collects estate-tax revenue $R_E$ from $\tau_e(b; \delta)$ on the bequest schedule (App. §D), where $\delta$ parameterizes the tax-rate shift applied to the SQ schedule under the reform. Under the status quo ($\delta = 0$ and *Verschonungsregeln* in place), aggregate revenue is $R_E^{\mathrm{SQ}}$ — the model-computed analog of the published €13.3 bn for 2024 — and flows entirely to general government. Under the counterfactual reform, aggregate revenue decomposes:

$$T_e \;=\; T_e^{\mathrm{SQ}} \;+\; T_e^{\mathrm{grant}}$$

where $T_e^{\mathrm{SQ}}$ is held constant at the model-calibrated SQ stationary revenue level across all counterfactual scenarios — a comparability constraint ensuring that the general-government outflow does not vary across scenarios — and $T_e^{\mathrm{grant}}$ funds the grant. Aggregate consistency:

$$T_e^{\mathrm{grant}} \;=\; G \cdot N_{h_0}$$

for the *Grunderbe* experiment, with analogous identities for the other transfer variants. The tax-rate shift $\delta$ is solved for in an outer loop (App. §G) such that aggregate counterfactual revenue $R_E^{\mathrm{reform}}(\delta)$ delivers $T_e^{\mathrm{grant}}$ equal to the target. The grant size $G$ is therefore an experiment *input* (with two values tested: €20k and €200k); $\delta$ is the outer-loop output. Experiment 4 is an exception: it fixes $\delta = 0$ but closes the *Verschonungsregeln*, and reports the implied $G$ as endogenous output.

**Cross-sectional density.** The joint density $g(a, z_p, z_\varepsilon, h, t)$ evolves under

$$\partial_t g = -\partial_a[s\, g] + \tfrac{1}{2}\partial_{aa}(\sigma_r^2 a^2 g) + \mathcal{L}_{z_p}^\top g + \mathcal{L}_{z_\varepsilon}^\top g - \partial_h g - m(h)\, g + \mathcal{B}[g; \tau_e, G(g)]$$

where $\mathcal{B}$ depends on $g$ through mortality outflux (the inherited wealth) and is parameterized by the contemporaneous policy instruments $(\tau_e, G)$. In each counterfactual experiment, $G$ is the target input and $\delta$ (the rate-shift parameter inside $\tau_e$) is found in an outer loop to deliver aggregate consistency. Stationary equilibrium is the joint solution $(V, g^*, \delta^*)$ to the coupled HJB–KF system at the $\delta^*$ that satisfies $T_e^{\mathrm{grant}}(\delta^*) = G \cdot N_{h_0}$.

---

## 3. Calibration

The model is calibrated to Germany at annual frequency.

**Preferences and life cycle.** Discount rate $\rho = 0.04$. Relative risk aversion $\gamma = 2$ (baseline), with sensitivity over $\gamma \in [1.5, 3]$. Households enter at $h = 0$ (biological 22), retirement at $h_\mathrm{ret} = 43$ (biological 65), terminal age $H = 78$ (biological 100). Age-specific mortality $m(h)$ from Destatis life tables, gender-pooled (cleaned CSV in `data/derived/mortality_hazard.csv`).

**Labor-productivity process — FSS 2010 level estimates.** The two stochastic components are calibrated to Fuchs-Schündeln, Krueger & Sommer (2010), using their preferred *level* moment-based estimates (which they argue are more plausible than the first-difference estimates; see FSS pp. 42–44). Specifically:

- *Persistent component innovation variance:* $\sigma_\eta^2 = 0.016$ per year (FSS Figure 19 bottom panel, household earnings, level estimation, averaged over time).
- *Transitory component variance:* $\sigma_\varepsilon^2 = 0.18$ (FSS Figure 19 bottom panel, household earnings, level estimation, averaged over time — the series ranges 0.15–0.22 over 1984–2002; we take the midpoint).
- *Life-cycle profile $\psi(h)$:* quartic in age, calibrated to the FSS auxiliary regression of log earnings on a quartic in age, evaluated at biological ages 22–65 and rescaled so that $\overline{\psi(h)} = 0$ over the working life.

For embedding in the continuous-time HACT model, we approximate FSS's discrete-time random-walk-plus-iid specification with two Ornstein–Uhlenbeck processes:

- Persistent: $dz_p = -\theta_p z_p\, dt + \sqrt{2\theta_p\,\sigma_{z_p}^2}\, dW_p$ with $\theta_p = 0.02$ (corresponding to annual persistence $\rho_p \approx 0.98$, an approximation of FSS's unit root chosen for stationarity; consistent with Bayer & Juessen 2009 who directly estimate $\rho \approx 0.9$ for German wages). Stationary variance $\sigma_{z_p}^2 = \sigma_\eta^2 / (1 - \rho_p^2) \approx 0.40$.
- Transitory: $dz_\varepsilon = -\theta_\varepsilon z_\varepsilon\, dt + \sqrt{2\theta_\varepsilon\,\sigma_{z_\varepsilon}^2}\, dW_\varepsilon$ with $\theta_\varepsilon = 3$ (annual persistence $\rho_\varepsilon \approx 0.05$, effectively iid per year). Stationary variance $\sigma_{z_\varepsilon}^2 = \sigma_\varepsilon^2 = 0.18$.

This is a deliberate departure from FSS's exact specification, with the trade-off documented in App. §A.1.4. FSS themselves note (p. 43) that "the simple permanent vs. transitory shock decomposition employed here is insufficient to describe wage and earnings shocks faced by German households," so we treat their estimates as defensible primary anchors but not exact truth. Refinement via De Nardi, Fella & Paz-Pardo (2020)-style nonlinear earnings dynamics is left as a follow-up.

The full process is discretized as a $7 \times 5 = 35$-state continuous-time Markov chain (7 persistent states via Rouwenhorst, 5 transitory states via Tauchen on iid Gaussian), with the two components independent. Conversion to the continuous-time generator via matrix logarithm; details in App. §A.1.

**Return process.** Baseline $r_0 = 0.02$. Wealth-dependent return: $\kappa = 0.02$ above the 95th percentile of the wealth distribution, based on Fagereng et al. (2020) adapted to German institutional features. Idiosyncratic return volatility $\sigma_r = 0.15$. Wage growth absent (long-run stationary equilibrium).

**Bequest motive.** Warm-glow strength $\phi$ calibrated by indirect inference to match the aggregate bequest-to-wealth ratio observed in *Erbschaftsteuerstatistik* (approximately 1% annually) and the wealth share of households over 65 in SOEP and PHF (approximately 40%). Threshold $\underline{a}$ calibrated to the empirical share of decedents leaving zero (or below-exemption) bequests in *Erbschaftsteuerstatistik*. The two parameters $(\phi, \underline{a})$ are jointly identified by these three moments via Nelder–Mead minimization of an identity-weighted moment-deviation objective; the outer loop is App. §F.

**Inheritance kernel.** Fertility distribution by completed cohort fertility from Destatis (TFR 1.35, 2024 vintage). Equal sharing across surviving offspring as benchmark; share of childless decedents (approximately 20%) routed to the general grant pool, augmenting $T_e^{\mathrm{grant}}$ (and therefore reducing the required tax-rate shift $\delta$ to hit a given target $G$). Heir's initial $z_p$ drawn from the persistent component's stationary distribution; heir's initial $z_\varepsilon$ drawn from the transitory stationary distribution. Both are independent of parental terminal states in the benchmark (intergenerational income transmission shut down); a sensitivity with persistence $\pi_0(z_p^\mathrm{heir} \mid z_p^\mathrm{parent})$ is in App. §C.4.

**Status-quo policy.** The current *Erbschaftsteuergesetz*: *Steuerklasse* I exemption €400,000 per child per parent, rate schedule as published (App. §D.1). *Verschonungsregeln* 85% exemption for business wealth, step up to 100% for estates above €26m (approximation in App. §D.2). Business-wealth share $\beta(a)$ from *Erbschaftsteuerstatistik* tabulations. **No grant in the status quo** ($G_\mathrm{SQ} = 0$); estate-tax revenue accrues to the general fund (modeled as outside the household budget set, since it does not feed back into any household-level instrument).

**Validation targets.** The status-quo stationary distribution is validated against (i) Albers-Bartels-Schularick (2024) top wealth shares: top 1% ≈ 26%, top 10% ≈ 60%, bottom 50% ≈ 2% (2021 figures); (ii) SOEP- and PHF-implied wealth-by-age profile; (iii) aggregate bequest-to-wealth ratio from *Erbschaftsteuerstatistik* (≈1% lower bound, since most below-exemption transfers are not tabulated); (iv) aggregate household net wealth-to-GDP ratio ≈ 4.5×, combining Bundesbank Financial Accounts (financial side) with Destatis Vermögensbilanzen (non-financial side), excluding statutory pension wealth; (v) **total *Erbschaftsteuer* revenue under SQ schedule ≈ €13.3 bn for 2024 (Destatis PD25_320) — this is the level pinned as $T_e^{\mathrm{SQ}}$ in the aggregate-consistency constraint of §2.** Validation-target tolerances in App. §J.

> **Calibration-failure flag.** If $\phi$ has to take an implausibly large value to hit both the W/Y target (4.5×) and the bequest-to-wealth target jointly, the return process or income process needs to absorb more of the variance. In that case the calibration is rerun with $\sigma_r$ as a third free parameter; if even that fails, the calibration tension is reported transparently. A related diagnostic: the top-1% share target of 26% is sensitive to $\kappa$, $\sigma_r$, and $\phi$; a stationary distribution that hits W/Y but misses the top-1% share by more than 3 percentage points indicates that the heterogeneous-return calibration needs revisiting.

---

## 4. Status-quo characterization

Three exercises before policy experiments.

**Decomposition of inequality.** What fraction of the top one percent's wealth share is attributable to (i) the wealth-dependent return premium, (ii) the bequest motive and intergenerational transmission, (iii) idiosyncratic income risk (persistent vs. transitory separately), (iv) the *Verschonungsregeln* exemption? Each margin shut down in turn, both in pure partial equilibrium (the additional revenue is not rebated — held by general government) and under the paper's aggregate-consistency convention (additional revenue funds the grant at the same per-capita level as Experiment 1a). Both numbers reported; the consistency-preserving variant is the headline.

**Life-cycle profiles.** Average wealth, consumption, and bequest probabilities by age under the calibrated stationary distribution, compared to SOEP and PHF empirical profiles. Particular attention to the wealth-by-age hump, which is shaped jointly by the life-cycle earnings profile $\psi(h)$, retirement, and the bequest motive.

**Intergenerational mobility.** Rank–rank correlation between parental and offspring wealth under the benchmark calibration (intergenerational income transmission shut down), compared to the sparse German evidence (Schnitzlein 2016 and follow-ups). The persistence sensitivity (App. §C.4) reports the additional channel.

---

## 5. Policy experiments

Four experiments. Each is a stationary-equilibrium comparison plus a transition path from the status-quo distribution to the new stationary state. Transition algorithm in App. §E. **All experiments respect the aggregate-consistency constraint of §2: $T_e^{\mathrm{SQ}}$ (the outflow to general government) is held constant at the model-calibrated SQ revenue level, and the additional revenue $T_e^{\mathrm{grant}}$ funds the policy instrument.** Experiments 1, 2, and 3 take the grant size as input; the tax-rate-shift parameter $\delta$ is found in an outer loop. Experiment 4 fixes $\delta = 0$ but closes the *Verschonungsregeln*, treating the resulting grant size as endogenous output.

**Experiment 1: the *Grunderbe* (two grant-size variants).** Structural reform: (i) closure of the *Verschonungsregeln* business-wealth exemption, (ii) reduction of the *Steuerklasse* I exemption to €200,000, (iii) tax-rate shift $\delta$ applied uniformly to the marginal-rate schedule of the SQ *Erbschaftsteuer*. Additional revenue $T_e^{\mathrm{grant}}$ paid as a one-time grant at age $h_0 = 0$ (biological 18). Two sub-experiments:

- **1a: $G = €20{,}000$** (the Bach 2021 / SPD–Greens headline). Outer loop on $\delta$ such that $T_e^{\mathrm{grant}}(\delta) = €20{,}000 \cdot N_{h_0}$.
- **1b: $G = €200{,}000$** (an order-of-magnitude larger; intentional upper-bound stress test). Outer loop on $\delta$ such that $T_e^{\mathrm{grant}}(\delta) = €200{,}000 \cdot N_{h_0}$.

For 1b, infeasibility is a possible (and instructive) outcome: if no $\delta \le \bar\delta$ (a reasonable upper bound on rate shifts) delivers the target without driving the total marginal rate above 100% somewhere, the experiment reports this as a finding about the revenue capacity of *Erbschaftsteuer* reform alone. For 1a, the expected $\delta$ is modest (the headline reform package — *Verschonungsregeln* closure plus a top-rate increase of perhaps 10–20 percentage points — should comfortably deliver €20k per cohort).

**Experiment 2: estate-tax reform with retirement transfer.** Same structural reform as Experiment 1a and the $\delta$ from 1a; the additional revenue $T_e^{\mathrm{grant}}$ is rebated as a one-time lump-sum transfer at age $h_\mathrm{ret} = 43$ instead, with transfer size $T_\mathrm{ret} = T_e^{\mathrm{grant}} / N_{h_\mathrm{ret}}$. Isolates the timing channel — same revenue, different lifecycle placement.

**Experiment 3: estate-tax reform with annual income transfer.** Same structural reform and $\delta$ from Experiment 1a; the additional revenue is paid as a flat annual income transfer to every working-age household, with transfer size $T_\mathrm{annual} = T_e^{\mathrm{grant}} / N_\mathrm{working}$. Isolates the lump-sum-vs-flow channel — same lifetime present value (approximately), different liquidity profile.

**Experiment 4: *Verschonungsregeln* closure only.** Close *Verschonungsregeln* but leave the rate schedule and exemption at SQ levels (i.e., $\delta = 0$). The additional revenue $T_e^{\mathrm{grant}} = R_E^{\mathrm{reform}} - T_e^{\mathrm{SQ}}$ is then endogenous and pays for an endogenous age-18 grant $G = T_e^{\mathrm{grant}} / N_{h_0}$. This is the "minimum-reform" experiment: it asks what grant size the *Verschonungsregeln* closure alone can fund, and provides a natural reference point for the political feasibility of Experiment 1 at each target grant size.

**Reporting.** For each experiment (and each sub-experiment of Experiment 1), four panels.

*Steady-state distribution.* Top 1%, top 5%, top 10%, bottom 50% shares; Gini; P99/P50 and P50/P10 ratios. Compared in tabular form across experiments and the status quo. The realized tax-rate shift $\delta$ and total revenue decomposition $(T_e^{\mathrm{SQ}}, T_e^{\mathrm{grant}})$ are reported as outputs alongside the distributional metrics.

*Distributional density plots.* Pre- and post-reform stationary $g(a)$ on a log scale, with policy-induced mass shifts highlighted.

*Transition dynamics.* Path of top shares and bottom-50% share over the 50 years following the reform, computed by the Achdou et al. (2022) transition algorithm: backward calendar-time HJB sweep from the new stationary value function, then forward calendar-time KF integration from the pre-reform $g^*_\mathrm{SQ}$. The tax-rate shift $\delta$ is held fixed at the new-stationary-equilibrium value along the entire path (the policy is announced once and remains in place); the realized $T_e^{\mathrm{SQ}}$ and $T_e^{\mathrm{grant}}$ vary along the path as the distribution evolves, with a small accounting adjustment per time step to maintain $T_e^{\mathrm{SQ}}$ at its reference level (any over- or under-shoot is folded into $T_e^{\mathrm{grant}}$, so the grant size $G_t$ fluctuates mildly along the path). App. §E details this.

*Welfare.* Consumption-equivalent welfare changes by initial wealth quintile and age. Under CRRA, $\lambda(a, z_p, z_\varepsilon, h) = (V_\mathrm{new}/V_\mathrm{SQ})^{1/(1-\gamma)} - 1$. Aggregate welfare decomposed into level, insurance, and intergenerational components following Floden (2001); formulas in App. §H.

---

## 6. Sensitivity analyses

Robustness over the calibration parameters whose values are least firmly established.

The bequest-motive strength $\phi$ over a range covering the alternative De Nardi (2004) calibrations. The wealth-return premium $\kappa$ over a range covering the lower and upper estimates from Fagereng et al. and reasonable German alternatives. Risk aversion $\gamma \in \{1.5, 2, 3\}$. The persistent-vs-transitory split: cut $\sigma_\eta^2$ in half (taking the lower end of FSS's wage estimates) and double-check that policy conclusions are robust. The retirement replacement rate $\mathrm{repl} \in \{0.43, 0.55, 0.65\}$ (OECD gross / OECD net / generous). The inheritance kernel (equal sharing versus single heir; ergodic versus persistent $z_p$ transmission). The *Steuerklasse* I exemption threshold $E_I \in \{€100\text{k}, €200\text{k}, €400\text{k}\}$ — at fixed grant target $G$, varying $E_I$ shifts which households contribute the additional revenue and therefore changes the required $\delta$; reported as a $\delta$-decomposition and a distributional-incidence comparison.

**Numerical robustness.** Three checks reported in a separate small table: (i) doubling the wealth-grid resolution $I$ changes top-1% shares by less than 0.2 pp; (ii) increasing $a_{\max}$ by 50% changes top-1% shares by less than 0.1 pp and probability mass in the top 1% of the grid stays below $10^{-5}$; (iii) halving the transition time step changes the top-1% path at $t = 10$ by less than 0.2 pp.

A single robustness table summarizes how the headline top-1% and bottom-50% share effects change across the economic specifications.

---

## 7. Expected mechanisms and findings

The paper should articulate, in advance of producing numbers, what mechanisms are expected to dominate and what surprises might arise.

**Expected dominant mechanisms.** The *Grunderbe*'s mechanical effect on the bottom of the distribution is large in the short run but dampened in the long run by household consumption of the grant. The headline finding for Experiment 1a (€20k) should be a substantial increase in the bottom-50% share and a moderate reduction in the top-1% share, achievable with a politically plausible $\delta$. The reform's contribution to reducing top concentration operates primarily through closure of the *Verschonungsregeln*, which mechanically taxes the wealth held in the form most concentrated at the top — Experiment 4 isolates this margin.

**Possible surprises.** First, the €200k stress test (Experiment 1b) may turn out infeasible — i.e., no $\delta$ within reasonable bounds delivers the required $T_e^{\mathrm{grant}}$. If so, the report quantifies *how* infeasible: what marginal-rate increase would be needed in principle, and what share of the bequest base would be taxed away. Second, heterogeneous returns may amplify or attenuate the grant's effect depending on whether recipients can access the higher-return wealth tier. Third, the transition dynamics may be slower than the policy debate assumes: bottom effects bind quickly (one generation), but top effects depend on cumulative estate-tax incidence over multiple generations.

**Robustness of conclusions.** The qualitative finding — that the *Grunderbe* plus *Verschonungsregeln* closure substantially reduces concentration but is fundamentally revenue-constrained — should be robust across calibrations. Quantitative magnitudes are sensitive to return-heterogeneity parameters and bequest-motive strength; this sensitivity is reported transparently.

---

## 8. Implementation and timeline

**Software stack.** MATLAB (R2021a or later). Base architecture from *schaab-lab/SparseEcon*, specifically `use_cases/08` (life-cycle / OLG) and `use_cases/05` (transition dynamics). Moll's reference codes at *benjaminmoll.com/codes/* and the *HACT_Numerical_Appendix.pdf* for technique. Code released open-source at publication.

**Module-to-codebase mapping.** Complete table in App. §I. The genuinely new pieces:

- The inheritance kernel $\mathcal{B}$ with target $G$ and outer-loop tax-rate shift (App. §C, §G);
- The German *Erbschaftsteuer* tax-and-exemption function with *Verschonungsregeln* (App. §D);
- The two-component income process (persistent + transitory) with continuous-time generator (App. §A.1);
- The indirect-inference wrapper for $(\phi, \underline{a})$ (App. §F);
- The Floden welfare decomposition (App. §H).

Everything else is adaptation of existing modules.

**Timeline (nine months).**

*Months 1–2.* SparseEcon setup, `use_cases/08` baseline run. Implement wealth-dependent return $r(a)$ and idiosyncratic return diffusion. Add warm-glow terminal condition $W(a)$. Validate against Achdou et al. (2022) and De Nardi (2004) benchmarks (App. §J tests 1–2).

*Months 3–4.* Implement two-component income process (App. §A.1). Implement the inheritance kernel $\mathcal{B}$ and the *Erbschaftsteuer* tax function. Verify mass conservation (test 3) and aggregate revenue against *Erbschaftsteuerstatistik* (test 5). Calibrate $(\phi, \underline{a})$ via indirect inference. Validate the status-quo stationary equilibrium.

*Months 4–5.* Status-quo characterization (Section 4): inequality decomposition, life-cycle profiles, intergenerational mobility.

*Months 5–6.* Four policy experiments. Steady-state distributional metrics.

*Month 7.* Transition dynamics for each experiment. Sensitivity analyses including numerical-robustness checks (tests 6, 7, 8).

*Months 8–9.* Write up. First-draft circulation to Bartels and other DIW colleagues, revision, submission to DIW Discussion Paper series.

---

## 9. Outputs and venues

**Primary output.** DIW Discussion Paper, eventually journal article. Target venues: *Journal of Public Economics*, *European Economic Review*, *Quantitative Economics*, *Review of Economic Dynamics*.

**Secondary outputs.** Open-source MATLAB code on GitHub with reproducible calibration and experiment scripts, shipping with the App. §J unit tests so any subsequent user can verify the build. A short *DIW Wochenbericht* policy brief summarizing the quantitative implications for the German debate.

**Strategic value for the Walter-Benjamin proposal.** The paper does five things: it builds the HJB–KF–transition-dynamics infrastructure WP2 will reuse; produces a citable working paper in the proposal bibliography by submission time; demonstrates end-to-end implementation capability; establishes a working relationship with the host (Bartels) on a topic of mutual interest; speaks directly to the German *Grunderbe* policy debate. The paper has independent academic value regardless of the proposal outcome.

---

## 10. Technical appendix

The companion document *Technical Appendix for Implementation* contains the full implementation-ready specification: state-space discretization and boundary conditions (§A), HJB–KF stationary algorithm (§B), inheritance kernel as one explicit equation (§C), *Erbschaftsteuer* tax function with rate-shift parameter $\delta$ (§D), transition algorithm (§E), indirect-inference loop for $(\phi, \underline a)$ (§F), aggregate-consistency closure with outer-loop $\delta$ for the target $T_e^{\mathrm{grant}}$ (§G), Floden welfare decomposition (§H), module-to-codebase mapping (§I), unit-test plan (§J).

---

## References

Achdou, Y., Han, J., Lasry, J.-M., Lions, P.-L., & Moll, B. (2022). Income and Wealth Distribution in Macroeconomics: A Continuous-Time Approach. *Review of Economic Studies* 89(1), 45–86.

Albers, T. N. H., Bartels, C., & Schularick, M. (2024). Wealth and Its Distribution in Germany, 1895–2021. DIW Discussion Paper 2105.

Andreoni, J. (1990). Impure Altruism and Donations to Public Goods: A Theory of Warm-Glow Giving. *Economic Journal* 100(401), 464–477.

Bach, S. (2021). Grunderbe und Vermögensteuern können die Vermögensungleichheit verringern. *DIW Wochenbericht* 88(50), 807–815.

Bayer, C., & Juessen, F. (2009). The Life Cycle and the Business Cycle of Wage Risk: A Cross-Country Comparison. Working Paper. *[Source for the German persistence estimate $\rho \approx 0.9$ that motivates the OU approximation of FSS's unit root.]*

De Nardi, M. (2004). Wealth Inequality and Intergenerational Links. *Review of Economic Studies* 71(3), 743–768.

De Nardi, M., Fella, G., & Paz-Pardo, G. (2020). Nonlinear Household Earnings Dynamics, Self-Insurance, and Welfare. *Journal of the European Economic Association* 18(2), 890–926. *[Possible follow-up extension for richer earnings dynamics.]*

Drechsel-Grau, M., Peichl, A., Schmieder, J. F., Schmid, K. D., Walz, H., & Wolter, S. (2022). Inequality and Income Dynamics in Germany. CESifo Working Paper 9605. *[Secondary source for descriptive trends on more recent data.]*

Fagereng, A., Guiso, L., Malacrino, D., & Pistaferri, L. (2020). Heterogeneity and Persistence in Returns to Wealth. *Econometrica* 88(1), 115–170.

Floden, M. (2001). The Effectiveness of Government Debt and Transfers as Insurance. *Journal of Monetary Economics* 48(1), 81–108.

Fuchs-Schündeln, N., Krueger, D., & Sommer, M. (2010). Inequality Trends for Germany in the Last Two Decades: A Tale of Two Countries. *Review of Economic Dynamics* 13(1), 103–132. *[Primary source for the parametric persistent-transitory earnings-process decomposition.]*

Schaab, A., & Zhang, A. T. (2021). Dynamic Programming in Continuous Time with Adaptive Sparse Grids. Working Paper. Code repository: `schaab-lab/SparseEcon`.

Storesletten, K., Telmer, C. I., & Yaron, A. (2004). Cyclical Dynamics in Idiosyncratic Labor Market Risk. *Journal of Political Economy* 112(3), 695–717.
