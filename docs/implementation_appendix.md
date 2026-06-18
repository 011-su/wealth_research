# Technical Appendix for Implementation

*Companion document to* **Estate-Tax-Financed Capital Grants and the Long-Run German Wealth Distribution: A Quantitative Analysis** *(preliminary working paper outline).*

This document is the implementation-ready specification for the MATLAB build. Each section is self-contained and intended to be picked up module-by-module against the corresponding source code from *schaab-lab/SparseEcon* (specifically `use_cases/08` for life-cycle / OLG, `use_cases/05` for transition dynamics) and the Moll reference codes at *benjaminmoll.com/codes/*.

Notation: $a$ wealth, $z_p$ persistent log-productivity, $z_\varepsilon$ transitory log-productivity, $h$ age, $t$ calendar time (transition paths only). $V(a, z_p, z_\varepsilon, h)$ value function, $g(a, z_p, z_\varepsilon, h, t)$ joint density, $g^*$ stationary. Indices $(i, j, k, \ell)$ run over discretization of $(a, z_p, z_\varepsilon, h)$.

---

## A. State-space discretization and boundary conditions

### A.1 Grids

**Wealth grid.** $I = 500$ points on $[a_{\min}, a_{\max}]$ with $a_{\min} = 0$ (no borrowing in benchmark; sensitivity with $a_{\min} = -0.3 \bar y$ is run) and $a_{\max} = 1{,}000 \bar y$ where $\bar y$ is mean labor income. Nonuniform geometric spacing:

$$a_i = a_{\min} + \frac{(a_{\max} - a_{\min})(e^{\lambda(i-1)/(I-1)} - 1)}{e^\lambda - 1}, \quad i = 1, \ldots, I,$$

with $\lambda = 5$ (densifies near the borrowing constraint, stretches for the Pareto tail).

### A.1.1 Persistent income component $z_p$

The persistent log-productivity follows a continuous-time Ornstein–Uhlenbeck process:

$$dz_{p,t} = -\theta_p\, z_{p,t}\, dt + \sigma_{p,\mathrm{diff}}\, dW_{p,t},$$

with $\theta_p = 0.02$ (annual persistence $\rho_p = e^{-\theta_p} \approx 0.98$) and $\sigma_{p,\mathrm{diff}}^2 = 2\theta_p \sigma_{z_p}^2$, where the stationary variance $\sigma_{z_p}^2$ is pinned by the FSS innovation variance $\sigma_\eta^2 = 0.016$:

$$\sigma_{z_p}^2 = \frac{\sigma_\eta^2}{1 - \rho_p^2} = \frac{0.016}{0.0396} \approx 0.404.$$

Discretize with $N_p = 7$ states via **Rouwenhorst** (better than Tauchen at high persistence). The output is a $7 \times 7$ annual transition matrix $P_p$ and the grid $\{z_{p,1}, \ldots, z_{p,7}\}$ centered on zero. Convert to a continuous-time intensity matrix via matrix logarithm:

$$\Lambda_p = \log_m(P_p) / \Delta t, \quad \Delta t = 1 \text{ year},$$

using MATLAB's `logm`. Verify (a) row sums equal zero to within $10^{-10}$ and (b) all off-diagonals are non-negative. If (b) fails (rare for Rouwenhorst), fall back to $\Lambda_p = (P_p - I)/\Delta t$ and document.

### A.1.2 Transitory income component $z_\varepsilon$

The transitory log-productivity follows a fast-mean-reverting OU process:

$$dz_{\varepsilon,t} = -\theta_\varepsilon\, z_{\varepsilon,t}\, dt + \sigma_{\varepsilon,\mathrm{diff}}\, dW_{\varepsilon,t},$$

with $\theta_\varepsilon = 3$ (annual persistence $\rho_\varepsilon = e^{-3} \approx 0.05$, essentially iid each year) and stationary variance $\sigma_{z_\varepsilon}^2 = \sigma_\varepsilon^2 = 0.18$ (FSS earnings level estimate).

Discretize with $N_\varepsilon = 5$ states via Tauchen on the stationary Gaussian. Generate the $5 \times 5$ annual transition matrix $P_\varepsilon$ (which will be close to a rank-1 matrix because of the fast mean reversion: all rows nearly equal to the stationary distribution). Convert to $\Lambda_\varepsilon$ via matrix log, applying the same validity check as §A.1.1 — but note that a near-rank-1 $P_\varepsilon$ frequently has *no* valid Markov generator (the embedding problem), and $\log_m(P_\varepsilon)$ will then carry negative off-diagonals. The fallback $\Lambda_\varepsilon = (P_\varepsilon - I)/\Delta t$ is the *expected* path for the transitory component, not a rare exception; verify non-negative off-diagonals and zero row sums to $10^{-10}$ either way, and document which branch was taken.

### A.1.3 Combined income state

The two components are independent. The combined state lives on $N_p \times N_\varepsilon = 35$ joint states. The combined generator is the Kronecker sum:

$$\Lambda_y = \Lambda_p \otimes I_{N_\varepsilon} + I_{N_p} \otimes \Lambda_\varepsilon$$

(both terms have the same sparsity pattern; the result is a sparse $35 \times 35$ generator with row sums zero and non-negative off-diagonals).

### A.1.4 Why this approximation of FSS

FSS specify the persistent component as a random walk ($\rho = 1$) and the transitory as iid. The continuous-time approximation differs on two fronts:

- *Persistent:* OU with $\rho_p = 0.98$ instead of unit-root random walk. Rationale: a true random walk has no stationary distribution, which is incompatible with the standard HACT framework. Bayer & Juessen (2009) directly estimate $\rho \approx 0.9$ for German wages, suggesting the unit-root assumption is itself an approximation. Our $\rho_p = 0.98$ is between Bayer-Juessen and FSS, conservative on the side of FSS. Stationary variance is finite (≈ 0.40), close to the variance the random walk would accumulate over ~25 years (0.40 ≈ 25 × σ²_η = 25 × 0.016).
- *Transitory:* OU with $\rho_\varepsilon = 0.05$ instead of pure iid. Rationale: continuous time doesn't have a natural notion of "iid per year"; the fastest mean reversion that's still numerically stable is $\theta_\varepsilon \sim 3$, giving annual persistence near 0.05 — close enough to iid that consumption smoothing handles it the same way.

FSS themselves note (p. 43) that the simple persistent/transitory decomposition is incomplete. Refinement to nonlinear earnings dynamics (De Nardi, Fella & Paz-Pardo 2020) is a future extension.

### A.1.5 Life-cycle profile $\psi(h)$

Quartic in age, calibrated to FSS's auxiliary regression of log earnings on age polynomial. Functional form:

$$\psi(h) = \psi_0 + \psi_1 h + \psi_2 h^2 + \psi_3 h^3 + \psi_4 h^4$$

where $h$ is model age. The coefficients $\{\psi_0, \ldots, \psi_4\}$ are extracted from the FSS regression tables, evaluated over biological ages 22–65, with $\psi_0$ rescaled so that $\overline{\psi(h)} = 0$ over the working life (this gives $\mathbb{E}[y] \approx 1$ after Jensen, modulo the small correction from the quartic shape interacting with the constant variance).

If FSS's exact coefficients are not directly transcribable from the published tables, fit the quartic by regressing FSS Figure 2's average-wage-by-age profile (or equivalently the cross-sectional age-earnings curve from SOEP if accessible) on the quartic. Document the source in the calibration scripts.

### A.1.6 Earnings level

Combining all pieces:

$$y(z_p, z_\varepsilon, h) = \exp\!\big(\psi(h) + z_p + z_\varepsilon - \tfrac{1}{2}(\sigma_{z_p}^2 + \sigma_{z_\varepsilon}^2)\big)$$

with the Jensen term ensuring $\mathbb{E}[y \mid h] = \exp(\psi(h))$. The grand mean $\overline{\mathbb{E}[y \mid h]}$ over the working life is normalized to 1 by the $\psi_0$ shift in §A.1.5, so the wage $w$ is literally mean earnings.

**Units bridge to EUR.** Model income and wealth are denominated in units of mean annual gross earnings ($w = 1$). The *Erbschaftsteuer* exemptions and brackets (§D) and the grant $G$ are, by contrast, EUR-denominated, so they are only well-defined relative to a single scale factor $\bar y_{\mathrm{EUR}}$ = mean annual gross earnings in EUR (German full-time equivalent, $\approx €45{,}000$; pin to the calibration vintage). A model wealth level $a$ is $a \cdot \bar y_{\mathrm{EUR}}$ euros — e.g. the €400,000 *Steuerklasse* I exemption is $E_I \approx 8.9$ in model units, and $a_{\max} = 1000\,\bar y \approx €45$m. Store $\bar y_{\mathrm{EUR}}$ in the `par` struct and convert every EUR input to model units once at construction; never mix the two in the operators.

### A.1.7 Retirement

At $h \geq h_\mathrm{ret} = 43$, the productivity state collapses to a single deterministic value $y_\mathrm{ret} = \mathrm{repl} = 0.55$ (OECD net replacement rate, Germany). Implementation: the joint income generator $\Lambda_y$ is replaced by an absorbing structure at $h_\mathrm{ret}$ — all 35 working-age states transition deterministically to a single retirement state at the moment of retirement, after which there is no income stochasticity.

In the discrete state representation, the retirement state can be treated as an additional 36th state with absorbing dynamics, accessed only at $h \geq h_\mathrm{ret}$. This avoids redundant generator computation in the retirement region.

### A.1.8 Age grid

$N_h = 79$ points, $h = 0, 1, \ldots, 78$. Age advances at unit rate. Working life $h = 0, \ldots, 42$; retirement $h = 43, \ldots, 78$. Mortality hazard $m(h)$ from `data/derived/mortality_hazard.csv`.

### A.1.9 Total state size

$I \times N_p \times N_\varepsilon \times N_h = 500 \times 7 \times 5 \times 79 = 1{,}382{,}500$ during working life (plus a smaller absorbing retirement region). All operators are sparse; the full HJB linear system has of the order $10^7$ nonzeros. *This is larger than the earlier draft's $5 \times 10^5$; the transitory component costs a factor of ~5.*

**Budget for memory, not just time.** Each implicit step is a direct sparse factorization of a system this size, and the LU fill-in — not the nonzero count of $A$ — drives the footprint. On the prototype at $N \approx 1.7 \times 10^5$ (one-eighth of this) a single stationary solve already took ~1–2 minutes and exhausted a 16 GB machine, with repeated out-of-memory crashes during the multi-solve policy experiments. At $N \approx 1.4 \times 10^6$ plan for a workstation with $\ge 32$–64 GB, run the grid-robustness checks (test 6) at reduced $I$ first, and consider state ordering ($a$ slowest, to cut the aging/kernel bandwidth) or switching the KF solve to `bicgstab` before committing to the full dense-tensor grid. "Tens of seconds" is realistic only with ample RAM and an ordered, possibly iterative, solve.

### A.2 Boundary conditions

**HJB.**

- *Borrowing constraint at $a = a_{\min}$.* State-constraint BC (Achdou et al. 2022, App. A): force forward drift at $i = 1$ to $\max(s_{1,...}^F, 0)$; equivalently $V'(a_{\min}, \cdot) = u'(w y)$ at the corner if the constraint binds. SparseEcon `use_cases/08` already implements this.
- *Upper wealth at $a = a_{\max}$.* Reflecting (zero-flux): $s_{I,...}^F = 0$. Verify after computation that less than $10^{-5}$ of probability mass sits in the top 1% of the grid.
- *Terminal age at $h = H$.* $V(\cdot, H) = W(a) = \phi(a-\underline a)^{1-\gamma}/(1-\gamma)$ for $a \geq \underline a$, with a tiny constant floor for $a < \underline a$.
- *Retirement transition at $h = h_\mathrm{ret}$.* Continuity: $V(a, z_p, z_\varepsilon, h_\mathrm{ret}^-) = V(a, y_\mathrm{ret}, h_\mathrm{ret})$ for all $(z_p, z_\varepsilon)$. The 35-dim working-age value collapses to a 1-dim retirement value.

**KF.**

- *Wealth boundaries.* Zero-flux at both ends; automatically enforced by upwind construction.
- *Age $h = 0$.* Source term from the inheritance kernel $\mathcal{B}$ (see §C).
- *Retirement transition.* Mass at $h_\mathrm{ret}^-$ from all 35 working-age states aggregates into the retirement state at $h_\mathrm{ret}$.

**Mass-conservation invariant.** Total population stationary: $\int \mathcal{B}\, da\, dz_p\, dz_\varepsilon = \int m(h) g\, da\, dz_p\, dz_\varepsilon\, dh$. Verify; see §J test 3.

---

## B. HJB–KF stationary algorithm

### B.1 Implicit upwind HJB sweep over age

For the life-cycle model, age is a deterministic backward-progressing state. SparseEcon `use_cases/08` solves the HJB backward in age from $V(\cdot, H) = W$:

**Algorithm B.1 (HJB backward sweep).**

```
Input: terminal V(:,:,:,H) = W(a)   (size I × N_p × N_eps)
For h = H-1, H-2, ..., h_ret:
    Retirement region: single income state, deterministic y_ret.
    Solve scalar HJB (Algorithm B.1a) with implicit step.
For h = h_ret - 1, ..., 0:
    Working-life region: 35 income states.
    Initialize V^(0) = V(:,:,:,h+1) reshaped from retirement boundary
    For n = 0, 1, 2, ...:
        Compute upwind dV_F, dV_B over a.
        Apply state-constraint BC at i=1 and zero-flux at i=I.
        c_F = (dV_F)^(-1/gamma); c_B = (dV_B)^(-1/gamma); c_0 = ra + wy.
        s_F = r(a)a + wy - c_F; s_B = ...; choose upwind per (i,j,k).
        Assemble sparse generator A^(n) of size (I·N_p·N_eps)^2:
            - drift block (tridiagonal in i, diagonal in (j,k))
            - diffusion block (tridiagonal in i) for sigma_r^2 a^2
            - persistent income block (Kronecker: Lambda_p ⊗ I_{N_eps})
            - transitory income block (Kronecker: I_{N_p} ⊗ Lambda_eps)
        Implicit step (Delta = 1000):
            (1/Delta + rho + m(h)) V^(n+1) - A^(n) V^(n+1) =
                u(c) + (1/Delta) V^(n) + m(h) W + (1/Delta_h) V(:,:,:,h+1)
        Solve linear system.
        If ||V^(n+1) - V^(n)||_inf < tol_V: break
    V(:,:,:,h) = V^(n+1)
    Store c*, s*, A.
```

Notes:
- The implicit relaxation parameter $\Delta \in [1000, 10000]$ is false-time, not real-time.
- Mortality term $-m(h)(V - W)$ split: $-m(h) V$ in the implicit operator, $+m(h) W$ on the RHS.
- Tolerance $\text{tol}_V = 10^{-8}$, typically 50–200 inner iterations per age with naïve relaxation.

**Howard's algorithm acceleration.** Replace the inner relaxation loop with policy iteration: alternate between (a) policy evaluation — one sparse solve of $(\rho + m(h)) V - A^c V = u(c) + m(h) W$ with the current consumption policy $c$ fixed, and (b) policy improvement — one update $c \leftarrow (\partial_a V)^{-1/\gamma}$ from the FOC. Convergence typically in 5–15 iterations instead of 50–200, a 5–10× speedup. Critical for the transition algorithm where the HJB is resolved at every calendar time step. Lift the wrapper from Moll's `HACT_Additional_Codes.pdf` reference code if SparseEcon `use_cases/08` doesn't already implement it.

### B.2 KF stationary solver

After the HJB sweep, the optimal drift $s^*$ and consumption $c^*$ are known. The KF stationary density solves

$$\mathcal{A}^*[g^*] + \mathcal{B}[g^*; \tau_e, G(g^*)] - m(h)\, g^* = 0, \quad \int g^* = 1,$$

with the adjoint generator $A^{*\top}$ obtained by transposing the assembled upwind operator from B.1 (dropping the consumption diagonal term).

### B.3 Coupled stationary solver with inheritance kernel and target grant

The inheritance kernel $\mathcal{B}$ takes the grant $G$ as input. The KF stationary density depends on $\mathcal{B}$, and the implied aggregate revenue $R_E^{\mathrm{reform}}$ depends on $g^*$ and the tax-rate shift $\delta$. The aggregate-consistency constraint requires $R_E^{\mathrm{reform}}(\delta) - T_e^{\mathrm{SQ}} = G \cdot N_{h_0}$. There are therefore two loops:

- **Inner KF loop** (this section): given $(\delta, G)$, solve for the stationary density $g^*$ under the contemporaneous kernel.
- **Outer δ loop** (§G): given target $G$, bisect on $\delta$ to satisfy aggregate consistency.

**Algorithm B.3 (KF stationary at given $(\delta, G)$).**

```
Inputs: tax-rate shift delta, target grant size G, structural reform flag
        (e.g., Verschonungsregeln_closed = true/false)
Initialize g^(0) = uniform on grid, normalized
For m = 0, 1, 2, ...:
    Compute B[g^(m); tau_e(delta), G] via §C
    Solve A*^T g^(m+1) - m(h) g^(m+1) + (ageing) + B^(m) at h=0 = 0
        with mass constraint sum(g) = 1
    If ||g^(m+1) - g^(m)||_1 < tol_g: break
```

Tolerance: $\text{tol}_g = 10^{-9}$. Typically 10–30 inner iterations. $G$ is held fixed throughout — it does not adjust to balance the budget here; that's the outer loop's job.

### B.4 Full stationary equilibrium

For Experiments 1, 2, 3, the equilibrium requires an outer loop on the tax-rate shift $\delta$ to satisfy aggregate consistency. For Experiment 4 (no rate change, Verschonungsregeln closed), $\delta = 0$ and the loop reduces to a single forward computation.

**Algorithm B.4 (full stationary equilibrium with target $G$).**

```
Inputs: target grant size G, structural reform parameters (E_I, Verschonung_closed)
Precompute: T_e^SQ_ref = revenue from SQ schedule applied to calibrated SQ g*
             (this is the model's analog of the published €13.3 bn; done once
              during calibration in §F)

Outer loop (bisection on delta):
    Initialize delta_lo = 0, delta_hi = delta_max
        (delta_max = 0.70, the hard cap at which the 30% top bracket reaches
         100%; see §G.7. Use the full [0, 0.70] bracket so the feasibility
         tests of §G.7 can be detected. A politically plausible sub-range is
         delta <~ 0.40, used only as a reference band in reporting, not as a
         search bound.)
    Repeat:
        delta = (delta_lo + delta_hi) / 2
        Run B.1 (HJB) under tax schedule tau_e(delta) with structural reform.
            If G enters HJB directly (Experiments 2, 3), use the current G.
            For Experiment 1, G enters HJB only via B, so the HJB doesn't
                actually depend on delta or G — solve once outside the loop.
        Run B.3 (KF stationary) at (delta, G).
        Compute R_E^total(delta) from §D.3 against g*(delta).
        Compute T_e^grant(delta) = R_E^total(delta) - T_e^SQ_ref.
        G_implied(delta) = T_e^grant(delta) / N_h0(delta).
        If G_implied(delta) < G_target:
            delta_lo = delta
        Else:
            delta_hi = delta
        If |G_implied(delta) - G_target| / G_target < tol_delta = 1e-3:
            break

For Experiment 4 (no rate change):
    delta = 0, Verschonung_closed = true
    Run B.1, B.3 once.
    G_endogenous = T_e^grant(0) / N_h0 — reported as output.
```

**Feasibility check.** For each target $G$, if $\delta_\mathrm{hi} = \delta_\mathrm{max}$ converges without ever delivering $G_\mathrm{implied} \ge G_\mathrm{target}$, the target is infeasible — report $\delta_\mathrm{max}$ as the rate shift that yields the largest grant achievable and quantify the shortfall. This is the expected outcome for the $G = €200{,}000$ stress test.

**A subtle point on coupling.** In Experiment 1 ($G$ paid at $h = 0$), $G$ enters the HJB only through the value of being born — i.e., through the initial distribution implied by $\mathcal{B}$. The HJB itself, conditional on $(a, z_p, z_\varepsilon, h)$, does not depend on $G$ directly. So a single HJB solve at the new tax schedule (which only depends on $\delta$) suffices; the outer δ-loop only re-runs B.3. In Experiments 2 and 3 (transfers at later ages or annually), the transfer enters the wealth drift and the HJB *does* depend on $G$ — those experiments require re-solving B.1 in each δ-iteration.

---

## C. Inheritance kernel $\mathcal{B}[g; \tau_e, G]$

### C.1 Continuous-form equation

$$
\bigl(\mathcal{B}[g; \tau_e, G]\bigr)(a_0, z_{p,0}, z_{\varepsilon,0}, 0)
= \int_0^H \int\!\!\int\!\!\int m(h)\, g(a, z_p, z_\varepsilon, h)\, \sum_{n \geq 0} f_n(n \mid h)\, K_n\, da\, dz_p\, dz_\varepsilon\, dh
$$

with $f_n(n \mid h)$ the offspring-count distribution conditional on parental age (Destatis), and the per-decedent kernel $K_n$:

- *For $n \geq 1$:*
$$K_n = \delta\!\left(a_0 - \frac{b(a; \tau_e)}{n} - G\right) \pi_p(z_{p,0}) \pi_\varepsilon(z_{\varepsilon,0}) \cdot n$$
where the multiplicative $n$ accounts for $n$ heirs receiving the same parcel, $b(a; \tau_e)$ is the after-tax bequest from §D, and $\pi_p, \pi_\varepsilon$ are the stationary distributions of the two income components (heirs draw $z_p$ and $z_\varepsilon$ independently from their ergodic distributions — intergenerational income transmission shut down in benchmark).

- *For $n = 0$ (childless decedent):* wealth $b(a; \tau_e)$ is added to the general grant pool (see §G).

Population stationarity: total mass injected by $\mathcal{B}$ at $h = 0$ equals total mortality outflux $\int m(h) g \, dh$.

### C.2 Discretization via the Young (2010) lottery

The Dirac mass $\delta(a_0 - b(a)/n - G)$ is resolved by the Young lottery. For each grid point $(a_i, z_{p,j}, z_{\varepsilon,k}, h_\ell)$ and each $n$:

1. Compute target $a^*_{i,n} = b(a_i; \tau_e)/n + G$.
2. Find indices $i', i'+1$ such that $a_{i'} \leq a^*_{i,n} \leq a_{i'+1}$.
3. Lottery weights $\omega_{i'} = (a_{i'+1} - a^*_{i,n})/(a_{i'+1} - a_{i'})$, $\omega_{i'+1} = 1 - \omega_{i'}$.
4. Contribution to $\mathcal{B}_{i'', j_0, k_0, 0}$:

$$\mathcal{B}_{i'', j_0, k_0, 0} \mathrel{+}= m(h_\ell)\, g_{i,j,k,\ell}\, f_n(n \mid h_\ell)\, \omega_{i''}\, \pi_{p,j_0}\, \pi_{\varepsilon,k_0}\, n$$

for $i'' \in \{i', i'+1\}$.

**Implementation note.** The Young (2010) lottery primitive — given a target value and a grid, produce two indices and two weights — appears throughout Moll's transition-dynamics code under names like `lottery` or `project_onto_grid`. Lift the primitive rather than rewriting; bugs in the lottery silently misallocate inherited wealth across the wealth grid and are hard to catch in inspection. The *application* to the inheritance kernel (the four nested loops, the heir-count multiplication, the income-state factorization) is genuinely new and must be written from scratch.

### C.3 Implementation as a sparse operator

Precompute $\mathcal{B}$ as a sparse linear operator $B$ of size $(I \cdot N_p \cdot N_\varepsilon) \times (I \cdot N_p \cdot N_\varepsilon \cdot N_h)$ that maps the flattened $g$ to the $h = 0$ source vector. Reconstruct once per outer iteration of B.3 (the operator only changes when $\tau_e$ or $G$ changes). Total nonzeros: $O(I \cdot N_p^2 \cdot N_\varepsilon^2 \cdot N_h \cdot n_{\max})$ with $n_{\max} = 4$ (maximum heir count). Roughly $10^7$ nonzeros — sparse but not tiny.

Note that since $\pi_p$ and $\pi_\varepsilon$ are *fixed* (the stationary distributions don't depend on parental state in the benchmark), the operator factors: $B = B_{wealth} \otimes \pi_p \otimes \pi_\varepsilon$. This factorization can reduce storage by a factor of $N_p N_\varepsilon = 35$. Implement the factored form; full assembly only for testing.

### C.4 Sensitivity: persistent income transmission

Benchmark: $\pi_p(z_{p,0}) = \pi_p^{\mathrm{erg}}$ (ergodic).

Sensitivity: $\pi_p(z_{p,0} \mid z_p^\mathrm{parent}) = P_p^k(z_{p,0} \mid z_p^\mathrm{parent})$, the $k$-step transition of the persistent process applied over a one-generational lag (e.g., $k = 30$ years). Calibrated to match an intergenerational income elasticity of 0.32 for Germany (Schnitzlein 2016). The transitory $z_\varepsilon$ remains ergodic in the sensitivity.

---

## D. *Erbschaftsteuer* tax function

### D.1 Rate schedule (current law, *Steuerklasse* I — children)

Taxable amount: $\tau_{\mathrm{base}}(a) = \max(0, (1 - \beta_{\mathrm{tax\text{-}base}}(a))\, a - E_I)$, where $\beta_{\mathrm{tax\text{-}base}}(a)$ is the fraction of the estate exempt under *Verschonungsregeln* (see D.2), and $E_I$ is the *Steuerklasse* I exemption ($€400{,}000$ under SQ; $€200{,}000$ under benchmark reform).

Marginal-rate schedule (*Steuerklasse* I), applied to $\tau_{\mathrm{base}}$. The SQ schedule is on the left; the reform introduces a uniform additive shift $\delta \in [0, \delta_\mathrm{max}]$ on every bracket, so the reform rate in bracket $\ell$ is $\min(\mathrm{rate}_\ell^\mathrm{SQ} + \delta, 1)$ (the cap at 100% prevents pathological cases):

| Taxable amount (€) | SQ marginal rate | Reform marginal rate |
|---|---|---|
| 0 – 75,000 | 7% | $\min(7\% + \delta, 100\%)$ |
| 75,001 – 300,000 | 11% | $\min(11\% + \delta, 100\%)$ |
| 300,001 – 600,000 | 15% | $\min(15\% + \delta, 100\%)$ |
| 600,001 – 6,000,000 | 19% | $\min(19\% + \delta, 100\%)$ |
| 6,000,001 – 13,000,000 | 23% | $\min(23\% + \delta, 100\%)$ |
| 13,000,001 – 26,000,000 | 27% | $\min(27\% + \delta, 100\%)$ |
| Above 26,000,000 | 30% | $\min(30\% + \delta, 100\%)$ |

Tax owed: $T_I(\tau_{\mathrm{base}}; \delta) = \sum_\ell \mathrm{rate}_\ell(\delta) \cdot \max(0, \min(\tau_{\mathrm{base}}, \mathrm{top}_\ell) - \mathrm{bottom}_\ell)$.

After-tax bequest: $b(a; \tau_e(\delta)) = a - T_I(\tau_{\mathrm{base}}(a); \delta)$.

Under SQ, $\delta = 0$ and $E_I = €400{,}000$ and *Verschonungsregeln* in place (per D.2). Under the structural reform, $\delta$ is the outer-loop variable in B.4, $E_I = €200{,}000$, and *Verschonungsregeln* are closed (Experiments 1, 2, 3); Experiment 4 fixes $\delta = 0$ and $E_I = €400{,}000$ but closes *Verschonungsregeln*.

*Steuerklassen* II and III implemented analogously but not used in the benchmark.

### D.2 *Verschonungsregeln* (business-wealth exemption)

Business-wealth share by estate-size bracket $\beta(a)$ taken from *Erbschaftsteuerstatistik* tabulations in `data/derived/erbsch_vermoegensart_shares.csv`. Parametric fit:

$$\beta(a) = \beta_\infty \cdot \frac{a^\eta}{a^\eta + a_\beta^\eta}$$

with starting values $\beta_\infty = 0.85$, $a_\beta = €5{,}000{,}000$, $\eta = 1.5$; refined against the data.

Under status quo: $\beta_{\mathrm{tax\text{-}base}}(a) = 0.15$ for $a < €26\text{m}$, $0.0$ for $a \geq €26\text{m}$ (the 85% / 100% step). Under Experiment 1/4 counterfactual: $\beta_{\mathrm{tax\text{-}base}}(a) \equiv 1$ for all $a$.

### D.3 Aggregate revenue

$$R_E(\delta) = \int\!\!\int\!\!\int\!\!\int m(h)\, g^*(a, z_p, z_\varepsilon, h)\, T_I\bigl(\tau_{\mathrm{base}}(a); \delta\bigr)\, da\, dz_p\, dz_\varepsilon\, dh$$

evaluated against the stationary $g^*$ implied by the model with the corresponding tax schedule. Define:

- $T_e^{\mathrm{SQ-ref}} \equiv R_E(\delta = 0)$ under the SQ schedule (no rate shift, *Verschonungsregeln* in place, $E_I = €400{,}000$), computed against the calibrated SQ stationary distribution. This is the model's analog of the published €13.3 bn. Computed once during calibration and stored.
- $R_E^{\mathrm{reform}}(\delta) \equiv R_E(\delta)$ under the reformed schedule (rate shift $\delta$, *Verschonungsregeln* closed, $E_I = €200{,}000$), computed against the counterfactual stationary distribution $g^*(\delta)$.
- $T_e^{\mathrm{grant}}(\delta) = R_E^{\mathrm{reform}}(\delta) - T_e^{\mathrm{SQ-ref}}$ — additional revenue funding the grant. Aggregate-consistency constraint requires $T_e^{\mathrm{grant}}(\delta) = G \cdot N_{h_0}(\delta)$.

Validation target (unit test 5): $T_e^{\mathrm{SQ-ref}}$ within 15% of €13.3 bn (Destatis 2024).

---

## E. Transition-dynamics algorithm

### E.1 Setup

Given pre-reform $(V^*_\mathrm{SQ}, g^*_\mathrm{SQ})$ and post-reform $(V^*_\mathrm{new}, g^*_\mathrm{new}, \delta^*)$ (the latter being the converged outer-loop value from B.4 at the target $G$), the transition is $\{V_t, g_t\}_{t=0}^T$ under the new tax policy $\tau_e(\delta^*)$, with:

- $V_T = V^*_\mathrm{new}$ (terminal),
- $g_0 = g^*_\mathrm{SQ}$ (initial),
- $\delta = \delta^*$ held fixed for all $t \in [0, T]$ — the policy is announced once at $t = 0$ and remains in place.
- The realized grant size $G_t$ along the path is determined by aggregate consistency at each $t$: $T_e^{\mathrm{grant}}(t) = R_E^{\mathrm{reform}}(\delta^*; g_t) - T_e^{\mathrm{SQ-ref}}$, $G_t = T_e^{\mathrm{grant}}(t)/N_{h_0}(g_t)$. In stationary equilibrium $G_t \to G_\mathrm{target}$; along the path $G_t$ fluctuates mildly as the distribution adjusts.

$T = 100$ years, $N_t = 200$ steps, $dt = 0.5$ year.

### E.2 Algorithm

**Algorithm E.1 (transition dynamics at fixed $\delta^*$, aggregate-consistency accounting per step).**

```
Inputs: V_new_star, g_SQ_star, delta_star (from B.4), T_e_SQ_ref (constant)

For q = 0, 1, 2, ... (outer iteration only for Experiments 2, 3):

    # Backward calendar-time HJB at fixed delta_star
    V[:,:,:,:,N_t] = V_new_star
    For t = N_t - 1, ..., 0:
        Run single implicit calendar-time step of B.1 with tax schedule
        tau_e(delta_star). For Experiments 2/3 use the path G_t computed
        below for the wealth-drift term; for Experiment 1 the HJB doesn't
        depend on G_t, only on delta_star.
        Store policies c*[:,:,:,:,t], s*[:,:,:,:,t], generator A[:,:,:,:,t].

    # Forward calendar-time KF with per-step aggregate-consistency accounting
    g[:,:,:,:,0] = g_SQ_star
    For t = 0, 1, ..., N_t - 1:
        R_E_total_t = revenue from g[:,:,:,:,t] under tau_e(delta_star) per §D.3
        T_e_grant_t = R_E_total_t - T_e_SQ_ref
        N_h0_t = ∫ g[:,:,:,:,t](:,:,:,h=0) da dz_p dz_eps
        G_t = T_e_grant_t / N_h0_t        # realized along the path
        Compute B[g[:,:,:,:,t]; tau_e(delta_star), G_t] via §C
        Forward step:
            g[:,:,:,:,t+1] = g[:,:,:,:,t]
                + dt * (A*^T[:,:,:,:,t] g - m(h) g + B_t + ageing)
        Renormalize sum(g) = 1.

    # Convergence check (only for Experiments 2, 3)
    If max_t change in V and g below tolerances: break
```

For Experiment 1, the loop converges in one pass because the HJB doesn't depend on $G_t$. For Experiments 2 and 3, typically 3–5 outer iterations.

**Note on the SQ-revenue constant.** $T_e^{\mathrm{SQ-ref}}$ is held fixed at its calibrated value throughout the transition. The path is therefore *not* in exact aggregate consistency at each $t$ if interpreted strictly: $T_e^{\mathrm{SQ}}$ (the share routed to general gov) might mechanically differ from $T_e^{\mathrm{SQ-ref}}$ when applied to a non-stationary distribution. The convention here is that the general-government allocation is held at the reference level by accounting fiat, and any mismatch is absorbed into $T_e^{\mathrm{grant}}$ (i.e., $G_t$). At the new stationary equilibrium this mismatch vanishes by construction.

### E.3 Reporting

For each transition: path of (top-1%, top-10%, bottom-50%, Gini, realized grant $G_t$) at annual frequency for $t = 0, \ldots, 50$. Plot overlay across the four experiments.

---

## F. Indirect-inference loop for $(\phi, \underline a)$

### F.1 Target moments

- $m_1$: bequest-to-wealth ratio. Target 0.01.
- $m_2$: wealth share of 65+. Target 0.40.
- $m_3$: share of decedents leaving below-exemption bequests. Target 0.30.

Model-implied moments computed from $g^*$ via straightforward aggregation.

### F.2 Objective and optimization

$$Q(\phi, \underline a) = \sum_{\ell=1}^3 \left(\frac{\hat m_\ell - m_\ell^\mathrm{target}}{m_\ell^\mathrm{target}}\right)^2$$

Minimize via `fminsearch` (Nelder–Mead). Each evaluation = one full call to Algorithm B.4. Typical: 30–60 inner solves; ~20–60 minutes wall time. Starting values $\phi = 10$, $\underline a = €100{,}000$. Stop when $Q < 10^{-3}$ or 100 evaluations.

### F.3 Failure handling

If $Q > 10^{-2}$ at the optimum, release $\sigma_r$ as third free parameter. If that fails, report calibration tension transparently per outline §3 calibration-failure flag.

---

## G. Aggregate-consistency closure with outer-loop $\delta$

This section specifies how the grant is funded under the aggregate-consistency convention introduced in outline §2. The estate-tax revenue under the reform decomposes as $T_e = T_e^{\mathrm{SQ}} + T_e^{\mathrm{grant}}$ where $T_e^{\mathrm{SQ}}$ is held constant at the SQ reference level for comparability across scenarios, and only $T_e^{\mathrm{grant}}$ funds the grant. For target grant size $G$, an outer loop on the tax-rate shift $\delta$ delivers aggregate consistency. **There is no income-tax instrument anywhere.**

### G.1 The SQ reference

$T_e^{\mathrm{SQ-ref}}$ is the model-computed aggregate estate-tax revenue under the SQ schedule (no rate shift, *Verschonungsregeln* in place, $E_I = €400{,}000$) applied to the calibrated SQ stationary distribution:

$$T_e^{\mathrm{SQ-ref}} \;=\; R_E(\delta = 0;\, g^*_\mathrm{SQ}) \;\approx\; €13.3 \text{ bn}$$

(Destatis target). Computed once at the end of the SQ calibration (§F) and stored as a constant. Validation: unit test 5.

### G.2 The aggregate-consistency constraint

For target grant $G$, the counterfactual stationary equilibrium satisfies:

$$\underbrace{R_E^{\mathrm{reform}}(\delta;\, g^*(\delta))}_{\text{total revenue under reform}} \;=\; \underbrace{T_e^{\mathrm{SQ-ref}}}_{\text{constant outflow to general gov}} \;+\; \underbrace{G \cdot N_{h_0}(\delta)}_{\text{grant pool}}$$

Rearranging: $T_e^{\mathrm{grant}}(\delta) \equiv R_E^{\mathrm{reform}}(\delta) - T_e^{\mathrm{SQ-ref}} = G \cdot N_{h_0}(\delta)$.

The outer loop in Algorithm B.4 solves for $\delta^*$ that satisfies this constraint at the target $G$.

### G.3 Childless decedents

The wealth of childless decedents (not claimed by heirs; share approximately 20%, per Destatis fertility data) is folded into the grant pool. Define

$$R_{\text{childless}}(\delta) \;=\; \int\!\!\int\!\!\int\!\!\int m(h)\, g^*\, f_0(h)\, b(a; \tau_e(\delta))\, da\, dz_p\, dz_\varepsilon\, dh$$

where $f_0(h)$ is the childlessness rate at parental age $h$. The grant identity becomes:

$$G \cdot N_{h_0}(\delta) \;=\; T_e^{\mathrm{grant}}(\delta) \;+\; R_{\text{childless}}(\delta)$$

(So the outer loop bisects $\delta$ to satisfy $T_e^{\mathrm{grant}}(\delta) + R_{\text{childless}}(\delta) = G \cdot N_{h_0}(\delta)$, with $T_e^{\mathrm{grant}}$ defined as above.)

### G.4 Stationary case

For each experiment with target $G$, Algorithm B.4 returns $\delta^*$ and the converged $g^*(\delta^*)$. Report:

- $\delta^*$ (the required rate shift)
- $T_e^{\mathrm{SQ-ref}}$ (held constant)
- $T_e^{\mathrm{grant}}(\delta^*) = G \cdot N_{h_0}(\delta^*) - R_{\text{childless}}(\delta^*)$
- Total revenue $R_E^{\mathrm{reform}}(\delta^*) = T_e^{\mathrm{SQ-ref}} + T_e^{\mathrm{grant}}(\delta^*)$

Status quo: $\delta = 0$, $G_\mathrm{SQ} = 0$, all revenue $T_e^{\mathrm{SQ-ref}}$ accrues to general government (modeled as outside the household budget set; does not feed back to households).

### G.5 Transitional case

Along the transition path, $\delta^*$ is held fixed at its stationary value. $T_e^{\mathrm{SQ-ref}}$ is also held fixed (by accounting convention; see §E.2). The realized $T_e^{\mathrm{grant}}$ and grant $G_t$ vary with the contemporaneous distribution but converge to their stationary values as $g_t \to g^*(\delta^*)$.

### G.6 Transfer variants for Experiments 2 and 3

Same $\delta^*$ as Experiment 1a (G = €20k) — i.e., the additional revenue pool $T_e^{\mathrm{grant}}$ is held at the level of Experiment 1a — but rebated differently:

**Experiment 2 (retirement transfer):**
$$T_\mathrm{ret} \;=\; \frac{T_e^{\mathrm{grant}} + R_{\text{childless}}}{N_{h_\mathrm{ret}}}$$
Paid as a one-shot wealth addition at $h = h_\mathrm{ret}$. Enters the HJB through the wealth-level jump $a^+ = a^- + T_\mathrm{ret}$ at the retirement boundary.

**Experiment 3 (annual income transfer):**
$$T_\mathrm{annual} \;=\; \frac{T_e^{\mathrm{grant}} + R_{\text{childless}}}{N_\mathrm{working}}$$
Paid annually to every working-age household. Enters the HJB through the wealth drift: $s = r(a)a + w y + T_\mathrm{annual} - c$ for $h < h_\mathrm{ret}$.

**Experiment 4 (Verschonungsregeln closure only):** fixed $\delta = 0$, $E_I = €400{,}000$, *Verschonungsregeln* closed. The grant is endogenous:
$$G \;=\; \frac{T_e^{\mathrm{grant}}(\delta = 0; \text{Verschonung closed}) + R_{\text{childless}}}{N_{h_0}}$$
Reported as model output rather than input.

### G.7 Feasibility of the $G = €200{,}000$ target

The bisection on $\delta$ may fail to converge if no $\delta \le \delta_\mathrm{max}$ delivers enough revenue. Specifically:

- $\delta_\mathrm{max}$ defined as the rate shift at which the highest bracket reaches 100% — for SQ top of 30%, $\delta_\mathrm{max} = 70\%$.
- If even at $\delta_\mathrm{max}$ the revenue $R_E^{\mathrm{reform}}(\delta_\mathrm{max}) - T_e^{\mathrm{SQ-ref}} < G \cdot N_{h_0}$, the target is infeasible.
- For $G = €200{,}000$ this is the expected outcome (the required aggregate revenue ~€150 bn against a bequest base of ~€400 bn/year would imply effective tax rates above feasibility).
- Report: maximal feasible $G$ at $\delta = \delta_\mathrm{max}$, plus the implied effective rate on top-bracket bequests. This is a substantive finding about the revenue-raising capacity of *Erbschaftsteuer* reform alone.

---

## H. Welfare decomposition

### H.1 Consumption-equivalent welfare change

Under CRRA:

$$\lambda(a, z_p, z_\varepsilon, h) = \left(\frac{V_\mathrm{new}(a, z_p, z_\varepsilon, h)}{V_\mathrm{SQ}(a, z_p, z_\varepsilon, h)}\right)^{1/(1-\gamma)} - 1$$

(care with signs: $V$ negative for $\gamma > 1$).

Report $\lambda$ averaged over $g^*_\mathrm{SQ}$ (utilitarian aggregate) and within initial wealth quintile × age bins.

### H.2 Floden (2001) decomposition

Aggregate $\bar\lambda$ decomposed:

$$\bar\lambda = \lambda^\mathrm{level} + \lambda^\mathrm{insurance} + \lambda^\mathrm{ig}$$

with:
- $\lambda^\mathrm{level} = \log(\bar c_\mathrm{new}/\bar c_\mathrm{SQ})$
- $\lambda^\mathrm{insurance} = \frac{1}{2}\gamma(\sigma^2_\mathrm{SQ} - \sigma^2_\mathrm{new})$ where $\sigma^2$ is cross-sectional log-consumption variance
- $\lambda^\mathrm{ig}$ = residual (intergenerational redistribution component)

**Implementation note.** The Floden decomposition appears as standard post-processing in several of Moll's published applications (e.g. the HACT JEEA paper's welfare computations). The formula is the same regardless of model — only the inputs ($V_\mathrm{new}$, $V_\mathrm{SQ}$, $g^*_\mathrm{SQ}$, $c^*_\mathrm{new}$, $c^*_\mathrm{SQ}$) change. Lift the skeleton from Moll's reference welfare code and plug in our state-space-specific aggregators.

---

## I. Module-to-codebase mapping

| Module | Source | Modification | Spec |
|---|---|---|---|
| Wealth grid | SparseEcon `lib/` | None | §A.1 |
| Rouwenhorst $z_p$ | SparseEcon `lib/` (write if absent) | New for HACT — Rouwenhorst not always in default util | §A.1.1 |
| Tauchen $z_\varepsilon$ | SparseEcon `lib/` | None | §A.1.2 |
| Matrix-log generator | New (one-liner with `logm`) | New | §A.1.1–2 |
| Kronecker-sum $\Lambda_y$ | New (`kron`) | New | §A.1.3 |
| Life-cycle profile $\psi(h)$ | New (FSS extraction) | New | §A.1.5 |
| Retirement collapse | New | New | §A.1.7, §A.2 |
| Mortality hazard | New (Destatis CSV) | None | §A.1.8 |
| HJB upwind FD | SparseEcon `use_cases/08` | Extend to 4D state $(a, z_p, z_\varepsilon, h)$ | §A.2, §B.1 |
| **HJB Howard's algorithm acceleration** | Moll `HACT_Additional_Codes` (if SparseEcon doesn't have it) | Wraps the inner iteration; 5–10× speedup, critical for transition algorithm | §B.1 |
| Sparse generator assembly | SparseEcon `lib/` | Extend with Kronecker income block | §B.1 |
| CRRA utility + FOC inversion | SparseEcon `lib/` (utility module) | None | §B.1 |
| Warm-glow $W$ | New (one line) | New | §A.2 |
| KF stationary (no kernel) | SparseEcon `lib/` | None | §B.2 |
| Young (2010) lottery primitive | Moll transition codes | None | §C.2 |
| **Inheritance kernel $\mathcal{B}$ (logic)** | **New (Young primitive lifted)** | **Build (§C); factored sparse form** | **§C** |
| ***Erbschaftsteuer* tax fn with $\delta$ shift** | **New** | **Build (§D); piecewise schedule, $\delta$ rate-shift parameter, *Verschonung* multiplier** | **§D** |
| Fertility distribution | New (Destatis CSV) | None | §C.2 |
| KF coupled (with kernel) | Extension of SparseEcon | Inner KF at given $(\delta, G)$ | §B.3 |
| **Outer-loop bisection on $\delta$ for target $G$** | **New** | **Aggregate-consistency closure: $R_E^{\mathrm{reform}}(\delta) - T_e^{\mathrm{SQ-ref}} = G \cdot N_{h_0}$** | **§B.4, §G** |
| **Aggregate-consistency revenue accounting** | **New** | **$T_e^{\mathrm{SQ-ref}}$ stored once; $T_e^{\mathrm{grant}}$ computed per equilibrium** | **§D.3, §G** |
| **Indirect-inference loop for $(\phi, \underline a)$** | **New** | **`fminsearch` wrapper; runs B.4 at each evaluation** | **§F** |
| Transition algorithm | SparseEcon `use_cases/05` | Adapt to 4D state; fixed $\delta^*$; per-step revenue accounting | §E |
| Distribution statistics (Gini, top shares, Lorenz, percentiles) | SparseEcon `lib/` + Moll Aiyagari examples | None — nonuniform-grid percentile computation is bug-prone, reuse tested code | §J reporting, all experiments |
| Welfare decomposition (Floden) | Moll Floden welfare example for skeleton | Adapt aggregators to our state space | §H |
| Plotting | New (MATLAB), templates from Moll figures where applicable | None | — |
| Unit-test harness | New | Per-test scripts | §J |

Bold rows are genuinely new code. Everything else adapts existing modules. Net assessment: the inheritance kernel logic, the *Erbschaftsteuer* tax function with the $\delta$-shift parameter, the outer-loop bisection on $\delta$ (with associated revenue accounting), and the indirect-inference wrapper are the pieces with no direct precedent in the cited codebases. Everything else is adaptation, lifting, or trivial.

---

## J. Unit-test plan

Each test is a separate MATLAB script in `tests/`. All tests must pass before calibration runs on real targets.

**Test 1 — Achdou et al. (2022) Huggett replication.** No life-cycle, no bequest, no return heterogeneity, no transitory income. Verify stationary asset distribution shape and equilibrium interest rate to 3 sig figs.

**Test 2 — De Nardi (2004) warm-glow replication.** Parameters from De Nardi Table 1 baseline. Top-1% wealth share within 0.5pp of her Table 4 (warm-glow case); W/Y within 5%.

**Test 3 — Mass conservation.** $|\sum g^* \Delta a\, \Delta z_p\, \Delta z_\varepsilon\, \Delta h - 1| < 10^{-10}$.

**Test 4 — Stationarity.** $\|A^{*\top} g^* + \mathcal{B}[g^*] - m \odot g^* + \mathrm{ageing}\|_\infty < 10^{-8}$.

**Test 5 — *Erbschaftsteuer* revenue.** Under status-quo schedule and calibrated $g^*$, $R_E$ within 15% of €13.3 bn (Destatis 2024).

**Test 6 — Wealth-grid robustness.** Doubling $I$ from 500 to 1000 changes the top-1% share in Experiment 1 by less than 0.2 pp.

**Test 7 — Upper-grid robustness.** Increasing $a_{\max}$ by 50% changes the top-1% share by less than 0.1 pp; mass in top 1% of grid stays below $10^{-5}$.

**Test 8 — Time-step robustness.** Halving $dt$ in transition algorithm changes top-1% path at $t = 10$ by less than 0.2 pp.

**Test 9 — Income process diagnostics.** Verify that the discretized two-component income process reproduces (within 2%):
- Persistent stationary variance $\sigma_{z_p}^2 \approx 0.40$
- Transitory stationary variance $\sigma_{z_\varepsilon}^2 \approx 0.18$
- $\mathbb{E}[y] = 1$ at every working age (after Jensen normalization)
- Annual first-order autocorrelation of $\log y$ matches the implied $\rho_p$ weighted by variance shares

**Test 10 — Sign and magnitude sanity (Experiment 1a, $G = €20{,}000$).** In Experiment 1a: (a) bottom-50% share strictly increases vs SQ; (b) top-1% share strictly decreases vs SQ; (c) outer loop converges to a finite $\delta^* \in [0, 0.4]$; (d) $T_e^{\mathrm{grant}}(\delta^*)$ within 1% of $G \cdot N_{h_0}$ (less the childless-decedent contribution).

**Test 11 — Aggregate-consistency invariant.** At the converged stationary equilibrium for any experiment, verify

$$\bigl| R_E^{\mathrm{reform}}(\delta^*) - T_e^{\mathrm{SQ-ref}} - G \cdot N_{h_0}(\delta^*) + R_{\text{childless}}(\delta^*) \bigr| / (G \cdot N_{h_0}) < 10^{-3}$$

i.e., the outer-loop convergence tolerance is met and the closure holds. For Experiment 4 (endogenous $G$), instead verify $G - (T_e^{\mathrm{grant}}(0) + R_{\text{childless}}(0))/N_{h_0}(0)$ is within $10^{-6}$ relative.

**Test 12 — Feasibility report for $G = €200{,}000$.** Run Experiment 1b. The expected outcome is that the bisection either (a) converges to $\delta^* \in (\delta_\mathrm{hi}^{\text{political}}, \delta_\mathrm{max})$ implying a politically extreme rate increase, or (b) fails to converge (infeasibility). The test is satisfied by either outcome — what's tested is that the algorithm correctly detects which.

Each test logged to `tests/results.log` with pass/fail and tolerance gap.

---

## Build order

Implement modules in this order, with each module gated by its unit test:

1. §A.1 wealth grid + §A.1.8 age grid + §A.1.5 life-cycle $\psi$ → no test
2. §A.1.1–3 income process (Rouwenhorst + Tauchen + Kronecker + matrix-log) → §J test 9
3. §B.1 HJB without $\mathcal{B}$ (4D state) → §J test 1
4. §B.2 KF without $\mathcal{B}$ → §J test 1
5. §A.2 terminal $W$ → §J test 2
6. §D tax function → no test (used in 7)
7. §C inheritance kernel + §G aggregate-consistency closure (outer-loop $\delta$) → §J tests 3, 4, 5, 11
8. §B.3 KF inner + §B.4 outer loop → §J tests 3, 4, 11
9. §F calibration loop → produces converged $(\phi, \underline a)$ and stores $T_e^{\mathrm{SQ-ref}}$
10. §E transition algorithm → §J tests 6, 7, 8
11. Experiments 1a (G=€20k), 1b (G=€200k), 2, 3, 4 → §J tests 10, 12
12. §H welfare decomposition → no test (sanity check on signs)

---

## Notes for implementation

- Read parameters from `data/derived/*.csv` and the FSS PDF in `data/papers/`. The FSS quartic coefficients for $\psi(h)$ should be extracted by a short script the first time the model is built; cache as `data/derived/psi_h_coeffs.csv`.
- All parameters in a single `par` struct passed by reference. Document units.
- Sparse `sparse(i,j,v,m,n)` constructors with preallocated triplets. No full matrices for $A$, $A^{*\top}$, $B$.
- MATLAB `\` for sparse solves; fall back to `bicgstab` only on profiling evidence.
- Save intermediate `g_star.mat`, `V_star.mat` after milestones (after B.1, B.3, B.4, calibration, each experiment).
- Each test script appends one line to `tests/results.log`: timestamp, test name, pass/fail, tolerance gap.
- No tab characters in `.m` files; 4-space indentation.

End of technical appendix.
