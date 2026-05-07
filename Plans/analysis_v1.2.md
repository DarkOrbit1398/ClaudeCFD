# Analysis v1.2: Path to 6th-order Compact Finite Difference 2nd-Derivative Operators

**Date:** 2026-05-07
**Sources:** Kim07, Kim24, Lele 1992, Jonathan Tyler Thesis, Chen 2023, Hixon CCD 2017, Chu & Fan CCD 1998, Deshpande 2021, Zhou & Zhang 2011, Brady & Livescu 2019, Dendro-GR 2022, Dendro-GR 2023, and supplementary web search.

This document is a **proposal for review**, not a decided plan. It supersedes analysis_v1.1.md and incorporates all corrections, clarifications, and new information from the review_followup.md exchange.

**Key corrections from v1.1:**

- §2.1 free-parameter argument was wrong; corrected below.
- Kim24 does not "abandon spectral error" for the interior; he retains Kim07's spectral optimization there. Our decision to use spectral optimization at boundaries is principled, not a deviation.
- Per-row $\hat{y}$ is not a deviation from Kim24; his paper is ambiguous on this point.
- Track B is not "invent wider stencils"; it is "diagnose why B-schemes underperform."
- Brady & Livescu 2019 does **not** use machine learning; see §4.3.

---

## 1. Background Recap (Corrected)

### 1.1 Kim07 (Interior Optimization)

Kim07 constructs optimized compact schemes for first derivatives using a spectral error functional. For an interior pentadiagonal-LHS / 7-point-RHS scheme with the standard antisymmetric parametrization

$$\beta f''_{i-2} + \alpha f'_{i-1} + f'_i + \alpha f'_{i+1} + \beta f''_{i+2} = \frac{1}{\Delta x}\!\left[a_1(f_{i+1}-f_{i-1}) + a_2(f_{i+2}-f_{i-2}) + a_3(f_{i+3}-f_{i-3})\right]$$

the scheme has five parameters $\{\alpha, \beta, a_1, a_2, a_3\}$. Two are consumed by the 4th- and 6th-order Taylor constraints (Kim07 Eqs. 2–3), leaving **three free parameters** for spectral optimization.

**Kim07's spectral error functional (Eq. 7):**

$$E = \int_0^r \Bigl[(1+\delta)\kappa\,D(\kappa) - N(\kappa)\Bigr]^2 \!\left(\frac{\kappa}{r}\right)^{\!n} d\kappa$$

where $N(\kappa) = 2(a_1\sin\kappa + a_2\sin 2\kappa + a_3\sin 3\kappa)$ and $D(\kappa) = 1 + 2\alpha\cos\kappa + 2\beta\cos 2\kappa$ are the RHS numerator and LHS denominator respectively, so that the modified wavenumber is $\bar{\kappa} = N/D$.

**Role of $n$ (power weight):** The factor $(\kappa/r)^n$ biases the integral toward higher wavenumbers, concentrating the optimization on short-wavelength resolution. With $n=10$, the high-$\kappa$ end of $[0, r]$ contributes much more to $E$ than the low-$\kappa$ end.

**Role of $\delta$ (small lift):** Kim replaces the target $\kappa$ by $(1+\delta)\kappa$ in the integrand. The paper describes this as "very effective if $\delta$ is adjusted properly." The value $\delta = -2.33\times10^{-4}$ is empirically determined via Mathematica iteration. It introduces a slight overall bias (a fractional upward/downward shift of the modified-wavenumber target) that minimally affects the resulting coefficients but is part of Kim07's fine-tuning.

**Relevant to our notebooks:** We tested the interior A4 scheme with and without $n$ and $\delta$ and found minimal difference in the resulting coefficients when compared against Kim07 Table 1. The parameters $n$ and $\delta$ were therefore discarded for all schemes (interior and boundary, first and second derivatives). This is an acceptable simplification.

**How the free parameters are determined:** Kim07 evaluates the stationary conditions $\partial E/\partial \alpha = 0$, $\partial E/\partial \beta = 0$, $\partial E/\partial a_1 = 0$ (i.e., choosing $\alpha$, $\beta$, $a_1$ as the optimization variables, with $a_2$, $a_3$ fixed by the Taylor constraints). The remaining two coefficients are then back-solved. The specific choice of which parameters to vary is empirical; Kim notes the results are insensitive to this choice (Table 2).

**Integration limit $r$:** Kim07 reports $r = 2.672 \approx 0.85\pi$, determined empirically by Mathematica iteration with no theoretical justification. This is an important point: the integration limit is not theoretically principled in Kim07; it is simply the value that produces the best spectral performance by trial. The role of $r$ is discussed in §2.2.

**First derivatives only.** Kim07 makes no claim about second derivatives.

### 1.2 Kim24 (Boundary Closure)

Kim24 retains Kim07's spectral-error interior and introduces a new boundary closure technique based on an extrapolating polynomial. The key clarification from v1.1: **Kim24 uses spectral optimization for the interior and physical-space $L_\infty$ error minimization for the boundary closure tuning only.** Kim24 Section 2.1 states explicitly that "the fourth-order interior scheme optimised by Kim [9] for maximised resolution performance is directly adopted," where Kim [9] is Kim07.

Our decision to use spectral optimization for both interior *and* boundary closures is intentional and theoretically motivated by the desire to optimize the operators independently of any specific test function. This is not a mistake.

**Boundary closure structure.** Kim24 introduces a degree-13 extrapolation polynomial $P(\bar{s})$ determined by 14 constraints (Eq. 10):

- 7 function-value matches: $P(j) = f_j$, $j = 0, \ldots, 6$
- 6 first-derivative matches: $\frac{dP}{d\bar{s}}(j) = f'_j$, $j = 0, \ldots, 5$
- 1 free constraint: $\frac{d^n P}{d\bar{s}^n}(\hat{y}) = 0$

The order $n$ in the free constraint must satisfy $n >$ (desired accuracy order) so that the constraint does not consume low-order Taylor content. For Kim24's 4th-order schemes, $n = 8$ is used. For our 6th-order schemes, we use $n = 9$ (with the observation that the polynomial coefficients become quartic in $\hat{y}$, balancing computation time against the number of local minima found).

**Per-row $\hat{y}$:** Kim24 Section 2.4 writes "once the free variables $\hat{y}$ are determined by the optimisation," using the plural. Upon close examination of Kim24's optimization tables and coefficient tables, the paper reports a single set of optimized boundary coefficients without specifying per-row values explicitly, leaving the question genuinely ambiguous. We should not get caught up on this detail: we have confirmed in our own work that using a separate $\hat{y}$ value for each boundary row works correctly for first derivatives, and there is no reason to change this approach.

Kim24 uses $\hat{y} \in (-6, 6)$ as the optimization range (Eq. 24), implemented via `scipy.optimize.differential_evolution`. Our current notebooks scan $\hat{y} \in [0, 5]$ in steps of 0.03; the range can be extended slightly into negative values (a few operators have local minima near $\hat{y} \sim -0.1$), but generally the minima are concentrated in the positive range.

**Non-Dirichlet weighted boundary probe.** Kim24 introduces a modified eigenvalue test via a matrix $\mathbf{Q}_r$ in which the first column is replaced by $(1-r)\mathbf{q}_0$ (Eqs. 16–18). This is a **boundary condition weighting parameter**, not a test function. With $r = 1$, the first column is zero (Dirichlet boundary: all weight on the enforced condition). With $r = 0.7$, the boundary node retains some influence of the interior solution, exposing eigenvalue instabilities that pure Dirichlet analysis hides. Kim24 Figure 2a shows that $r = 0.7$ moves eigenvalues further into the left half of the complex plane compared to $r = 1$, improving the stability margin. The value $r = 0.7$ was determined empirically by Kim24 "after some trial and error."

We are currently not using this probe (equivalent to $r = 1$). The implementation is straightforward: in `CFD Eigenvalue Stability.nb`, replace the first column of $\mathbf{Q}$ with $(1-r)\mathbf{q}_0$ before computing eigenvalues. A plan for when and how to implement this is in §5 (Track A).

**Stability constraints.** Kim24 enforces two stability constraints during optimization:

1. Eigenvalue constraint: all eigenvalues of $\mathbf{P}^{-1}\mathbf{Q}_r$ (with $r = 0.7$) must have $\text{Re}(\lambda) \geq 0$ (Kim24's sign convention; this corresponds to $\text{Re}(\lambda) \leq 0$ in our convention).
2. RK4 polynomial bound: $|R(-\text{CFL}\cdot\lambda_i)| < 1$ at CFL = 1.047.

We currently perform stability checking after-the-fact rather than during optimization. This is addressed in Track A.

**First derivatives only.** Kim24 says nothing about $f''$.

### 1.3 Our Extensions to Date

We have extended both the interior optimization (Kim07) and boundary closure (Kim24) frameworks to second derivatives by replacing the first-derivative polynomial constraint $\frac{dP}{d\bar{s}}(j) = f'_j$ with $\frac{d^2 P}{d\bar{s}^2}(j) = f''_j$ (confirmed in our notebooks; see §1.4 and §2.4 below for discussion of what this changes).

**Notebook state.** The current notebook coverage is approximately:

| Notebook | 1st deriv | 2nd deriv |
|---|---|---|
| `CFD Interior Optimization.nb` | T2, A4, B4, C4, A6, B6, C6, A8; hybrid schemes (Hybrid A42K1T, A41K2T) | A4, A6 (B4/C4/B6/C6 may be incomplete due to version overwrite) |
| `CFD Boundary Optimization.nb` | A4, B4, C4, A6, B6, C6, A8 | A4, B4, C4, A6, B6, C6 |
| `CFD Eigenvalue Stability.nb` | A4, B4, C4, A6, B6, C6, A8 | A4, B4, C4, A6, B6, C6 |

Known issues:

- 2nd-deriv P6(A6) `Output Coefficients` cell: `yHat1List[[5]]` indexes a length-4 list. Placeholder value left from editing; to be fixed.
- B4/C4/B6/C6 second-deriv work and hybrid interior schemes may have been lost to a Git version overwrite. The above table treats the current notebook as the authoritative state; once notebooks are restored, the table will be updated.
- S6 (heptadiagonal LHS, 1st deriv) is paused; suspected derivation typo. Not on the critical path.

**Per-row $\hat{y}$** confirmed working for first derivatives. We use separate values $\hat{y}_0, \hat{y}_1, \hat{y}_2$ for each boundary row, which provides more flexibility and is at worst equivalent to Kim24's approach.

**Hybrid interior schemes** (see §4) are a first-class option in the notebook; see Section `Spectral Schemes → Tyler Spectral Schemes` in `CFD Interior Optimization.nb`.

### 1.4 Our Spectral Error Functional

Because we discarded $n$ and $\delta$ from Kim07's functional, our actual spectral error integrands differ from Kim07 Eq. 7. The functional we minimize is:

**First-derivative interior** (minimized over $\kappa \in [0, r]$):

$$\Delta^2_{\!f'}(\kappa) = (\kappa\, D(\kappa) - N(\kappa))^2$$

with $N(\kappa) = 2(a_1\sin\kappa + a_2\sin 2\kappa + a_3\sin 3\kappa)$ and $D(\kappa) = 1 + 2\alpha\cos\kappa + 2\beta\cos 2\kappa$.

**Second-derivative interior** (minimized over $\kappa \in [0, r]$):

$$\Delta^2_{\!f''}(\kappa) = (\kappa^2 D(\kappa) - N(\kappa))^2$$

with $N(\kappa) = 4\!\left(a_1\sin^2\!\tfrac{\kappa}{2} + a_2\sin^2\!\kappa + a_3\sin^2\!\tfrac{3\kappa}{2}\right)$.

Note that $4\sin^2(m\kappa/2) = 2(1-\cos(m\kappa))$, so $N(\kappa)$ is a purely even function of $\kappa$ with no DC offset, as required for a symmetric second-derivative stencil.

**Boundary (first and second derivatives).** Because boundary spectral functions are complex-valued (capturing both dispersion and dissipation), the boundary integrand uses a complex norm. For second derivatives:

$$\Delta^2_{\text{bdy}}(\kappa) = \left(\kappa^2 - \left|\Re(\tilde{\kappa}(\kappa)) + i\,\Im(\tilde{\kappa}(\kappa))\right|\right)^2$$

The relevant code is in `CFD Boundary Optimization.nb` under any `Second Derivative → P6 (A6) Scheme → Optimize $\hat{y}$ → Node 2` section.

---

## 2. Diagnosis: Why Current 2nd-Derivative Operators Plateau

### 2.1 Free-Parameter Accounting (Corrected from v1.1)

The v1.1 claim that second-derivative A6 schemes have one fewer free parameter than first-derivative A6 schemes was **incorrect**, and is withdrawn. The correct count for both is:

| derivative | LHS params | RHS params | total | Taylor constraints | free for optimization |
|---|---|---|---|---|---|
| 1st-deriv A6 | $\alpha, \beta$ | $a_1, a_2, a_3$ | 5 | 2 | **3** |
| 2nd-deriv A6 | $\alpha, \beta$ | $a_1, a_2, a_3$ | 5 | 2 | **3** |

For first derivatives, the RHS is antisymmetric ($a_{-m} = -a_m$), so the zeroth-order consistency (stencil annihilates constants) is automatic. For second derivatives, the RHS is symmetric ($a_{-m} = a_m$) and $N(0) = 4\sum_m a_m \cdot 0 = 0$ identically — the zeroth-order constraint is also automatic. Both derivative orders have the same two independent Taylor constraints consuming two parameters.

The real source of the "Nate vs Tyler" error gap in the v1.1 smoking-gun table is not a lost free parameter but the **low-$\kappa$ vs total-$\kappa$ trade-off** described in §4. The table remains meaningful as evidence of this trade-off; it is just framed differently in v1.2.

### 2.2 The Integration Limit $r$ and the Role of the Spectral Integrand

#### 2.2a What $r$ controls

The integration limit $r$ is the **cutoff wavenumber** up to which the spectral error is minimized. Its effect is a direct trade-off:

- **Larger $r$**: the optimizer fits the modified wavenumber well over a wider range of $\kappa$, but the high-$\kappa$ tail ($\kappa$ near $\pi$) is forced into the fit, which tends to cause overshoot or oscillation in the modified wavenumber curve. Chen 2023 demonstrates explicitly that $r = \pi$ causes harmful overshoot for second-derivative schemes.
- **Smaller $r$**: the optimizer ignores high-$\kappa$ behavior, and the modified wavenumber curve naturally relaxes toward the standard (unoptimized) CFD scheme at high $\kappa$. Accuracy at low $\kappa$ is excellent but the scheme offers no improvement over standard CFD at high $\kappa$.

There is no universal principled rule for choosing $r$; it depends on what wavenumber range matters for the physical application. For numerical relativity, the relevant scales are set by the gravitational wave frequency relative to the grid spacing. As a practical guide: choose $r$ to be the highest wavenumber you need the scheme to resolve accurately. Frequencies above $r$ will be under-resolved regardless.

Chen 2023 recommends $r \in [0.5\pi, 0.75\pi]$ for second-derivative optimization, finding that $r \approx 0.85\pi$ (our current default) leads to overshoot in the modified wavenumber for second-derivative schemes. Our own tests (described below) support this.

#### 2.2b Our $r$-scan test results

You tested six values $r \in \{0.60\pi, 0.65\pi, 0.70\pi, 0.75\pi, 0.80\pi, 0.85\pi\}$ for A6 second-derivative operators, selecting approximately 3 operators per $r$ value (18 total). Results: 5 operators were stable in evolutions (2 at $r = 0.60\pi$, 1 at $r = 0.65\pi$, 2 at $r = 0.70\pi$); all stable operators converged at roughly **2nd order** rather than 6th order.

This is a significant finding that demands investigation. Two possible explanations:

1. **Boundary closure dominates the error.** If the boundary scheme is only 2nd-order accurate (or has a bug), the overall convergence order will be limited to 2nd order regardless of the interior scheme. With only 3 arbitrarily selected operators per $r$ value, we may have happened to select operators whose boundary coefficients are suboptimal or degenerate.

2. **Interior scheme quality degraded.** At smaller $r$, the interior modified wavenumber curve deviates significantly from $\kappa^2$ for mid-range $\kappa$. If the test function has significant energy at those wavenumbers, the convergence behavior could look lower-order.

The priority action is to test the full set of stable operators (not just 3 per $r$) and to isolate interior vs boundary contribution by testing with periodic boundary conditions (interior only) and then with boundaries.

#### 2.2c Scaling of the 2nd-derivative integrand

Our second-derivative integrand $(\kappa^2 D - N)^2$ is a degree-4 polynomial in $\kappa$ near $\kappa = 0$ (since both $\kappa^2 D$ and $N$ are $O(\kappa^2)$ at leading order), whereas the first-derivative integrand $(\kappa D - N)^2$ is degree-2 in $\kappa$. This means the 2nd-derivative integrand naturally downweights the low-$\kappa$ region relative to mid-$\kappa$. Even without the Kim07 $(\kappa/r)^n$ weight, the effective emphasis of the 2nd-derivative minimization is shifted toward higher $\kappa$ compared to the 1st-derivative case. This may explain why the 2nd-derivative optimization can show good low-$\kappa$ results on some metrics while struggling globally — the integrand itself does not strongly penalize low-$\kappa$ error.

#### 2.2d Role of $n$ and $\delta$ for 2nd derivatives

Since we have discarded $n$ and $\delta$, their direct role is moot for our notebooks. For completeness:

- Kim07's $(κ/r)^n$ with $n = 10$ was designed for first-derivative optimization where the target is linear in $\kappa$. For second derivatives the target is $\kappa^2$ and the integrand already carries an extra $\kappa^2$ factor. The combined effect of $(\kappa/r)^{10}$ on a $\kappa^4$-scaled integrand would dramatically over-emphasize high-$\kappa$, making the problem worse, not better.
- $\delta = -2.33\times10^{-4}$ is a 0.02% adjustment that has negligible effect on coefficients, as confirmed by your A4 test.

Conclusion: discarding $n$ and $\delta$ is correct for our purposes.

### 2.3 Boundary $\hat{y}$ Optimization

The current procedure — manual 1D scan of $\hat{y}$ for each boundary row, then combinatorial enumeration of stable combinations — is documented and functional. The main limitations are:

- **Range.** We scan $\hat{y} \in [0, 5]$ in steps of 0.03. Kim24 uses $(-6, 6)$. Extending slightly into negative values is warranted, particularly for operators where minima near $\hat{y} \sim -0.1$ have been observed.
- **Granularity.** The 0.03 step size has been empirically validated to produce smooth error curves with well-resolved minima. Finer resolution is unlikely to change results significantly.
- **Stability is post-hoc.** Accuracy is optimized first; stability is verified afterward by trimming. Accuracy and stability can conflict. A joint optimization using Mathematica's `NMinimize` with the stability constraint `Max[Re[Eigenvalues[M]]] ≤ 0` as an inequality is possible and is planned in Track A.

**The $n$ parameter in the zero-derivative constraint.** From your tests: for A6 schemes, choosing $n$ so that the polynomial coefficients are quartic in $\hat{y}$ is a good balance between runtime and number of minima found. Qualitatively, operators produced with different values of $n$ (e.g., $n \in \{8, 9, 10\}$) do not show meaningful differences in accuracy or stability — only in the number of candidate operators generated. This is consistent with Kim24's guidance that $n$ should exceed the desired order of accuracy; beyond that threshold, the specific value of $n$ is not critical.

### 2.4 Even/Odd Symmetry and 2nd-Derivative Spectral Fragility

First-derivative operators are **antisymmetric** under grid reflection: the stencil coefficients satisfy $Q_{i,j} = -Q_{i, 2i-j}$, and the interior modified wavenumber is purely real (no dissipation in the interior). Second-derivative operators are **symmetric**: $Q_{i,j} = Q_{i, 2i-j}$, and the interior modified wavenumber for a correct scheme is real and negative.

For boundary closures, both cases produce complex-valued modified wavenumbers (mixing dispersion and dissipation). The asymmetry between parity classes affects the spectral error function structure:

- **1st-deriv boundary:** the imaginary part of $\bar{\kappa}$ represents physical dissipation that can damp instabilities.
- **2nd-deriv boundary:** the sign of the real part of the modified wavenumber determines stability (must be non-positive; positive real part in $\tilde{\kappa}$ corresponds to a growing mode). The spectral error function for 2nd derivatives is therefore more sensitive to small coefficient perturbations that flip the sign of $\text{Re}(\tilde{\kappa})$ — this is the "spectral fragility" observed in practice.

Physically: for a second-derivative operator used in $u_{tt} = c^2 u_{xx}$, we need the spatial operator to behave like $-|k|^2$ modally. Any perturbation that makes $\text{Re}(\tilde{\kappa}) > 0$ at some wavenumber will produce an exponentially growing mode. The penalty landscape for avoiding this is stiffer than the corresponding landscape for 1st-derivative sign flips.

### 2.5 Non-Dirichlet Eigenvalue Probe

As clarified above (§1.2), the non-Dirichlet probe is a boundary condition weighting parameter in the eigenvalue matrix, not a test function. It is unrelated to the $L_\infty$ optimization in Kim24 and does not conflict with our preference for spectral optimization.

Our stability criterion for second derivatives uses eigenvalues $\mu$ of $\mathbf{P}^{-1}\mathbf{Q}$:

- Convention 1: $\mu = \mathbf{P}^{-1}\mathbf{Q}$ eigenvalues; stability requires $\text{Re}(\mu) \leq 0$ and $\text{Im}(\mu) = 0$.
- Convention 2: $\lambda = \sqrt{\mu}$; stability requires $\text{Re}(\lambda) = 0$ and $\text{Im}(\lambda) \geq 0$.

These are equivalent and the Mathematica code uses Convention 2 (the `eigs = (Eigenvalues[mat])^(1/2)` line).

Currently our eigenvalue analysis uses $r = 1$ (pure Dirichlet: first column of $\mathbf{Q}$ is zero). Switching to $r = 0.7$ should be straightforward; the full plan is in Track A below.

**Concerning the outlier eigenvalues.** We have observed operators with eigenvalue spectra that are correct in shape but feature outlier values with large imaginary part while still satisfying $\text{Re}(\lambda) \leq 0$. These have sometimes correlated with time-domain instabilities in wave-equation tests. The non-Dirichlet probe flattens the eigenvalue distribution toward the imaginary axis (see Kim24 Figure 2a); it is possible but unconfirmed that this diagnostic would expose such outliers as unstable. Testing with both $r = 1$ and $r = 0.7$ on operators that have produced outlier eigenvalues is the recommended first step.

---

## 3. Notebook State Summary

The v1.1 table is considered out of date pending restoration of overwritten notebook content (B4/C4/B6/C6 second derivatives and hybrid interior schemes). The table in §1.3 is the best current estimate.

Once notebooks are restored, an updated table should be produced. In the meantime, priority is on A6 second-derivative operators since this is the test case for all proposed improvements.

---

## 4. Hybrid Interior Schemes

### 4.1 Tyler's Approach vs Kim07

The Tyler thesis (Sections 3–4) develops spectral schemes with a fundamentally different optimization philosophy from Kim07:

- **Kim07** minimizes the integrated spectral error $E = \int_0^r (\ldots)^2 d\kappa$ via stationary conditions $\partial E / \partial \alpha = 0$, $\partial E / \partial \beta = 0$, $\partial E / \partial a_1 = 0$. This implicitly emphasizes the low-$\kappa$ region (since small $\kappa$ contributes most to the integral at our ranges).
- **Tyler** enforces exact spectral matching at selected wavenumbers: $\tilde{\omega}(\omega_i) = \omega_i$ for chosen $\omega_1, \omega_2, \omega_3$. These are interpolation-style constraints rather than integral-minimization conditions. Tyler distributes the matching across the full wavenumber range $[0, \pi]$, which tends to produce better high-$\kappa$ performance at the cost of losing the low-$\kappa$ guarantee.

Tyler's approach is described in the `Spectral Schemes → Tyler Spectral Schemes` section of `CFD Interior Optimization.nb`. Schemes labeled `P4 (A4)`, `P6 (A6)`, and `P6 (B6)` use Tyler-type constraints for all free parameters. The hybrids mix the two:

- **Hybrid A42K1T** (A4 scheme, 3 free params): 2 Kim-like + 1 Tyler-like constraint
- **Hybrid A41K2T** (A4 scheme, 3 free params): 1 Kim-like + 2 Tyler-like constraints

For A6 (3 free parameters), the analogous hybrids are `3K0T` (pure Kim), `2K1T`, `1K2T`, `0K3T` (pure Tyler).

### 4.2 The Trade-Off and What Hybrids Buy

The v1.1 "smoking gun" error gap — Nate's scheme winning low-$\kappa$ by 10× but losing total error by 6× against the Tyler baseline — is the expected consequence of using a pure Kim-like optimization. The hybrid schemes exist precisely to find a middle ground. This trade-off should be the central motivation for the interior-optimization work, replacing the (incorrect) "lost free parameter" framing.

Recommended next step: systematically vary the Kim/Tyler mix for A6 second-derivative interior schemes and plot error profiles (low-$\kappa$ error vs total spectral error). This characterizes the trade-off surface and identifies whether any hybrid combination closes the accuracy gap observed in evolutions.

### 4.3 Hybrid Weighting and the $\kappa^2$ Target (Clarification of D3)

The earlier v1.2 draft mentioned "the right hybrid weighting for $f''$" (item D3). To clarify: the hybrid approach (mixing Kim-like and Tyler-like constraints) applies identically to second derivatives — the same types of constraints work, applied to the second-derivative spectral error functional. The phrase "right weighting" was referring not to the Kim/Tyler constraint mix but to the $\kappa$-dependence of the integrand itself (item C3 / Zhou-Zhang weight). These are two separate questions:

1. **Constraint mix** (Kim/Tyler): Should be explored as hybrids for 2nd derivatives exactly as for 1st derivatives.
2. **Integrand form** (Zhou-Zhang weight): Could we replace $(\kappa^2 D - N)^2$ with $(\kappa^2 D - N)^2 \xi(\kappa)$ for some $\xi$ that shifts the optimization behavior? Addressed in §5 (Track A).

---

## 5. Proposed Plan

### Track A — Tighten What We Have (1–3 weeks, high-confidence gains)

**A1. Fix the `yHat1List[[5]]` index bug** in 2nd-deriv P6(A6). One-line fix; prerequisite for all other A-track work.

**A2. Implement the non-Dirichlet probe ($r = 0.7$).**

Step-by-step:
1. In `CFD Eigenvalue Stability.nb`, locate the section where the matrix $\mathbf{Q}$ is assembled.
2. Add a parameter `rND` (the non-Dirichlet weighting; set to 0.7 for Kim24's value, 1.0 for standard Dirichlet).
3. Replace the first column of $\mathbf{Q}$ with $(1-\text{rND}) \cdot \mathbf{q}_0$.
4. Re-run eigenvalue analysis for all current A6 2nd-deriv operators with `rND = 0.7`.
5. Compare plots: do outlier eigenvalues change character between `rND = 1` and `rND = 0.7`?
6. If the probe exposes previously "stable" operators as borderline, flag those for removal.

Expected difficulty: low. Expected benefit: may identify true instabilities masked by Dirichlet analysis.

**A3. Extend $\hat{y}$ scan range to $(-6, 6)$** for A6 2nd-deriv boundary optimization.

Step-by-step:
1. In `CFD Boundary Optimization.nb`, update the scan range from `Range[0, 5, 0.03]` to `Range[-6, 6, 0.03]`.
2. Re-run boundary optimization for A6 2nd-deriv, all three boundary rows.
3. Note any new local minima found in the negative range.
4. Enumerate the new full operator set and re-run stability analysis.

Expected difficulty: low. Expected benefit: may find stable operators missed by positive-only scan.

**A4. Investigate the convergence-order anomaly from the $r$-scan test.**

The finding of 2nd-order convergence for operators stable at $r \in \{0.60\pi, 0.65\pi, 0.70\pi\}$ needs diagnosis before any further $r$-scan work. Step-by-step:

1. Take one of the 5 stable operators and run a convergence test with *periodic* boundary conditions (isolating interior scheme quality).
2. Run the same test with physical boundaries using our current boundary closure.
3. Compare convergence order: if periodic gives 6th order but full gives 2nd order, the issue is in the boundary closure.
4. If periodic also gives 2nd order, the issue is in the interior coefficients at low $r$.
5. Plot the interior modified wavenumber curve for the selected operator; check whether it deviates from $\kappa^2$ significantly at mid-range $\kappa$.

This diagnosis should be done **before** running a larger operator sweep at reduced $r$.

**A5. Explore the Zhou-Zhang weight $\xi(\kappa)$ for the 2nd-derivative interior functional.**

Zhou & Zhang 2011 use $\xi(\kappa) = [1 + 2\alpha\cos\kappa + 2\beta\cos 2\kappa]^2$ as a weight in the second-derivative spectral error integral. This weight:
- Makes the integral analytically integrable (closed form, no numerical quadrature needed).
- Naturally emphasizes lower wavenumbers (where $\xi$ is smaller).
- **Does shift the location of optimal parameters** (confirmed by comparison of their Table 1 with unweighted coefficients).

Before adopting this weight, we should derive the stationary-point shift symbolically. For the first-derivative case, the shift is:

$$\frac{\partial}{\partial \alpha}\int_0^r (\kappa D - N)^2 \xi(\kappa)\,d\kappa = 0$$

vs

$$\frac{\partial}{\partial \alpha}\int_0^r (\kappa D - N)^2\,d\kappa = 0.$$

These produce different optimal $\alpha$ because $\xi$ is not constant over $[0, r]$. A symbolic computation in Mathematica will quantify whether the shift is significant. If the shift is small compared to the variation in optimal $\alpha$ across different $r$ values, then $\xi$ is safe to adopt for analytical convenience; otherwise, it changes the problem in a way that requires independent justification.

**A6. Implement `NMinimize` for joint $(\hat{y}_0, \hat{y}_1, \hat{y}_2)$ optimization.**

This replaces the current sequential 1D scans with a direct joint optimization. Step-by-step:
1. In `CFD Boundary Optimization.nb`, define a function `spectralError[y0_, y1_, y2_]` that computes the sum of boundary spectral errors for all three boundary rows given the three $\hat{y}$ values.
2. Call `NMinimize[{spectralError[y0, y1, y2], Max[Re[Eigenvalues[M[y0, y1, y2]]]] <= 0}, {{y0, -6, 6}, {y1, -6, 6}, {y2, -6, 6}}]`.
3. Compare the best operators found by this method against the best operators found by the current scan.

Note: the stability constraint makes this a constrained nonlinear program, which `NMinimize` handles natively. The computation may be slow (three-dimensional parameter space with eigenvalue evaluation), but a coarse initial run is feasible.

---

### Track B — Diagnose B-Scheme Failure (1–2 weeks)

Your description of the B-scheme failure: good spectral resolution, some stable operators, but convergence order lower than A/C schemes.

The B-schemes have a heptadiagonal RHS ($n_R = 4$, 9-point stencil), giving one additional free parameter over the A-schemes. Their poor convergence performance is unexpected if the interior is correct.

**Diagnosis steps:**

1. **Compare B6 interior modified wavenumber against A6.** Plot $\tilde{\kappa}(\kappa)/\kappa^2$ for the best B6 and A6 second-derivative interior operators on the same axes. If B6 deviates from the exact value ($= 1$) more than A6 at low or mid-$\kappa$, the interior coefficients have an accuracy issue.

2. **Check the Taylor constraint implementation in B6 interior.** The B-schemes add a parameter $a_4$ to the RHS. Verify that the Taylor constraints for 6th-order accuracy are correct for a 9-point RHS. A typo in the constraint equations would reduce the effective order without being obvious from spectral plots.

3. **Run isolated convergence tests (interior only, periodic BCs)** for B6 second-derivative operators. This isolates the interior accuracy from boundary effects.

4. **Compare B6 boundary closure coefficients** against what would be expected from the A6 boundary analysis. If the boundary notebooks for B6 use the wrong stencil parameters, the coefficients will be inconsistent with the interior.

5. **If all of the above pass**, the failure may be in how ghost zones are filled during evolution tests. Confirm that the B-scheme's wider RHS stencil is compatible with the ghost zone width used in the 1D evolution notebook.

---

### Track C — Combined Compact Difference (CCD) Approach (4–8 weeks)

#### C1. Why CCD

The key accuracy numbers from Chu & Fan 1998 (Table 1):

| scheme | $f'$ leading error | $f''$ leading error |
|---|---|---|
| 6th-order explicit central | $\sim 36/7!$ | $\sim 72/8!$ |
| Lele 6th tridiagonal compact | $\sim 4/7!$ | $\sim 17/8!$ |
| Lele 6th pentadiagonal compact | $\sim 5/7!$ | $\sim 28/8!$ |
| CCD (6th order) | $\sim 0.87/7!$ | $\sim 2/8!$ |

CCD reduces $f''$ leading error by **~8.4×** relative to the 6th-order tridiagonal compact scheme and **~13.8×** relative to the pentadiagonal compact scheme, on a **3-point stencil**. This is the most significant accuracy improvement available in the literature for second derivatives without going to radically wider stencils.

The mechanism is the joint treatment of $f'$ and $f''$ as unknowns at every node, coupled via a quintic Hermite polynomial (Chu & Fan Eqs. 2.7–2.8). Both derivatives are computed simultaneously in a single solve.

#### C2. Computational cost

**Linear system structure:** The CCD system is **twin-tridiagonal** — a block-tridiagonal system with 2×2 blocks (one row for $f'$, one row for $f''$ at each node). The solve is $O(N)$ via block forward/backward elimination (Chu & Fan Appendix 3), the same asymptotic cost as a standard scalar tridiagonal solve.

**Memory:** Storing both $f'$ and $f''$ at each node doubles the memory per field relative to standard compact FD. For Dendro-GR's 24 BSSN fields, this means storing 48 derivative arrays instead of 24. This is non-trivial but manageable.

**MPI cost:** In Dendro-GR, each block is 13×13×13 with $k = 3$ ghost zones per side (interior: 7×7×7). A CCD solve in one dimension operates on a 1D strip of 13 points (block size) at a time, with ghost zone data already communicated. The tridiagonal solve structure requires sweep from left to right and right to left within each strip, with no additional inter-block communication beyond what standard compact FD already requires. No additional ghost zones are needed compared to the current 6th-order explicit FD if the CCD stencil remains 3-point. This is a significant advantage over wider-stencil approaches.

**Summary:** CCD's cost overhead is the 2× memory per field and the twin-tridiagonal solve (approximately 2× the flop count of a single tridiagonal solve). Given the ~8–14× accuracy improvement for $f''$, this is likely worth the cost. A definitive decision on whether to pursue CCD should come after Track A results (to know whether we can achieve adequate $f''$ accuracy within the current pentadiagonal framework) and after completing the CCD derivation for our boundary case.

#### C3. Hixon 2017 Optimized CCD

Hixon 2017 extends CCD with two free parameters $\sigma_1$, $\sigma_2$ that modify the 4th-order error terms in the $f'$ and $f''$ equations respectively (Hixon Eqs. 5–6). He optimizes these over a fine $(1001 \times 1001)$ grid search, testing three objectives: minimize wave propagation error ($\alpha$) only, minimize dissipation error ($\beta$) only, and a blended criterion.

His **b5 boundary treatment** (the only one that achieves full design order in practice) uses five explicit stencils at each boundary, with the fifth stencil being a 9-point central difference that matches the first **eight** terms of the CCD Taylor expansion (Hixon Section IV.A). The b5 treatment is the correct template for CCD boundary implementation.

#### C4. CCD implementation plan

This is genuinely new derivation work. Step-by-step:

1. **CCD interior optimization (new notebook `CCD Interior Optimization.nb`):**
   a. Start from Chu & Fan Eqs. 2.7–2.8 (the coupled $f'$/$f''$ stencil).
   b. Introduce free parameters $\sigma_1$, $\sigma_2$ as in Hixon Eqs. 5–6.
   c. Derive the joint modified wavenumber for the coupled system.
   d. Define a spectral error functional targeting the $f''$ modified wavenumber (analogous to our current 2nd-deriv interior functional).
   e. Minimize with respect to $\sigma_1$, $\sigma_2$ and compare the resulting $f''$ performance against our best pentadiagonal results.

2. **CCD boundary closure (new notebook `CCD Boundary Optimization.nb`):**
   a. Extend Kim24's extrapolation polynomial to the **coupled** ($f'$, $f''$) system.
   b. The constraint structure changes: instead of one first-derivative condition per node, we impose *both* $\frac{dP}{d\bar{s}}(j) = f'_j$ and $\frac{d^2P}{d\bar{s}^2}(j) = f''_j$ at each node. This doubles the number of derivative conditions, requiring a higher-degree polynomial or a modified constraint count.
   c. Determine the polynomial degree needed for 6th-order CCD boundary closure (likely degree 13 is sufficient given the additional information available from $f''$).
   d. Free parameters $\hat{y}$ per boundary row as before; optimize as in Track A.

3. **CCD eigenvalue stability (update `CFD Eigenvalue Stability.nb`):**
   a. The CCD eigenvalue matrix is $2N \times 2N$ (two equations per node).
   b. The stability criterion for the joint operator is: eigenvalues of the $2N \times 2N$ system must satisfy the appropriate criterion (for the $f''$ subsystem, $\text{Re}(\lambda) = 0$ and $\text{Im}(\lambda) \geq 0$ in the square-root convention).
   c. Implement non-Dirichlet probe ($r = 0.7$) from the start.

4. **Validation:** Test on smooth functions (periodic BCs first, then with boundaries) and compare convergence order and error magnitude against best pentadiagonal A6 2nd-deriv results.

---

## 6. Ideas Considered and Not Recommended — Extended Justifications

### 6.1 Apply $D^1$ twice to get $D^2$

**Verdict: Skip.** We have direct empirical evidence that applying a first derivative twice gives worse accuracy than a dedicated second-derivative operator, and this is consistent with the analytical picture: applying $D^1$ twice introduces a composed modified wavenumber $\bar{\kappa}^2$ instead of the target $\kappa^2$, and the composition amplifies any spectral error in $D^1$ quadratically. The leading-error term for composed $D^1 \circ D^1$ is worse than a dedicated symmetric $D^2$ stencil of the same formal order.

### 6.2 Spotz-Carey HOC (High-Order Compact)

**Verdict: Skip.** Spotz & Carey 1996 achieve higher effective order by substituting the PDE (e.g., $u_{tt} = c^2 u_{xx}$) into the truncation error, canceling lower-order terms. This only works for a specific PDE. A generic second-derivative operator for use across BSSN, Maxwell, NSM, and Teukolsky cannot exploit PDE-specific cancellations. It is not a generic $D^2$ operator.

### 6.3 Ciment-Leventhal $(I + \delta^2/12)^{-1}\delta^2$

**Verdict: Skip.** This operator is algebraically equivalent to the standard 4th-order compact second-derivative scheme, but its derivation depends on the wave-equation structure ($u_{tt} = c^2 u_{xx}$). It is not useful as a standalone general-purpose $D^2$ operator for arbitrary fields.

### 6.4 Lee-Liu-Sun CCD2 (Mixed Derivatives, 2014)

**Verdict: Defer.** Their CCD extension handles mixed spatial derivatives ($\partial^2/\partial x \partial y$), which arise in full 3D BSSN but are not the current focus. If we implement CCD for pure second derivatives (Track C) and it performs well, the mixed-derivative extension is a natural follow-on. For now, the 1D second-derivative problem is not yet solved.

### 6.5 Central Compact Schemes (Liu/Shu 2013) and UTEP CHVM (2022)

**Verdict: Defer.** These schemes double the per-point degrees of freedom (storing values at both cell-center and cell-node points). Applied to BSSN's ~24 evolved fields, this doubles the memory footprint *on top of* the CCD 2× overhead, giving a ~4× total memory increase. The accuracy gains are real but the cost is high. If CCD alone provides adequate $f''$ accuracy, there is no reason to incur the additional cost of staggered/combined schemes.

### 6.6 Splitting Algorithms for High-Bandwidth Solves (Zhang 2016)

**Verdict: Not an accuracy tool.** Zhang 2016 addresses the communication cost of solving high-bandwidth (wider-than-tridiagonal) implicit systems on distributed hardware by splitting them into smaller local systems. This is a computational cost reduction technique. It would be relevant if we move to heptadiagonal LHS schemes (S6 track, currently deferred), but does not improve accuracy.

### 6.7 Zhou-Zhang Prefactorization (Zhou & Zhang 2011)

**Verdict: Worth knowing, not a priority.** The prefactorization decomposes an implicit compact scheme into forward and backward explicit biased stencils, eliminating the band-matrix solve at the cost of two explicit sweeps. The accuracy of the resulting scheme matches the original compact scheme's spectral resolution. This is a communication-efficiency win for production deployment (relevant for Dendro-GR's distributed architecture) but provides no accuracy improvement over standard compact FD. Worth revisiting at the production-code stage.

### 6.8 Deshpande KKT Framework (Deshpande 2021)

**Verdict: Worth prototyping as a tooling improvement.** Deshpande 2021 reformulates spectral-error optimization as a **quadratic programming problem with linear equality constraints**, solvable analytically via Lagrange multipliers (KKT conditions). For standard interior compact schemes on periodic domains, this gives closed-form optimal coefficients without iterative optimization. Key limitation: the KKT matrix can become rank-deficient at some stencil sizes, and the approach becomes more complex for non-periodic boundaries.

**Proposed prototype:** Apply the Deshpande framework to the A6 second-derivative interior scheme (a simple case with 3 free parameters) and compare the analytical result to our current `NMinimize` output. If the results match and the computation is faster, adopt it for interior optimization.

### 6.9 Brady & Livescu 2019

**Clarification from v1.1:** Brady & Livescu 2019 does **not** use machine learning. Their method is a direct nonlinear optimization that tests candidate boundary scheme parameters by running them on the compressible Euler equations and measuring stability of the solution. The optimization is a form of gradient ascent on a smoothness objective (Algorithm 1 in their paper). This is analogous to testing stability in physics rather than just checking eigenvalues.

Key features of their approach:
- Addresses 1st derivatives only (4th, 6th, 8th order).
- Constructs **conservative** boundary schemes (discrete conservation of globally conserved quantities).
- Enforces stability via direct nonlinear testing, not eigenvalue analysis.
- Produced operators that outperform several literature alternatives on nonlinear hyperbolic tests.

This is a useful reference for understanding what "stable boundary closure" means in a nonlinear setting. Their conservativity requirement may be relevant to our boundary closures in BSSN, where constraint damping involves conserved quantities. We should note their stability criterion (monotonicity of density in Euler; eigenvalue spectrum analysis) as a benchmark when we assess our own operators.

---

## 7. Additional Considerations

### 7.1 Dendro-GR Context

From the Dendro-GR papers (2022, 2023): Dendro-GR evolves BSSN with **24 fields** (1 lapse, 3 shift, 6 metric, 1 K, 10 trace-free extrinsic curvature, 3 conformal connections). The 2023 paper notes support for both BSSN and Z4c/CCZ4. Blocks are $13 \times 13 \times 13$ (confirmed), with $k = 3$ ghost zones per direction for 6th-order FD, leaving a $7 \times 7 \times 7$ interior.

**Constraints on CFD stencil size:**

- Maximum stencil radius: 3 grid points per direction (limited by $k = 3$ ghost zones).
- This is compatible with our current pentadiagonal-LHS / 7-point-RHS schemes (stencil radius 3).
- CCD (3-point stencil, radius 1) fits within this constraint easily.
- If we wanted to increase $k$, the papers indicate that the "unzip" (octant-to-patch) operation cost scales with ghost zone width, so increasing $k$ beyond 3 has a measurable computational cost.

**Implicit solve overhead in Dendro-GR:**

The papers do not discuss implicit solves (all time integration is explicit RK4). Adding compact FD would require 1D tridiagonal solves along each direction, which are local to each block. Within a block of 13 grid points per direction, a tridiagonal solve is trivially fast. The primary overhead is that the solve cannot be parallelized trivially across the direction of solve, but within a block this is not an issue. Inter-block CCD coupling would require communication of $f'$ and $f''$ values at block boundaries — the same ghost zone infrastructure already in place.

**Where CFD provides the most benefit:** The RHS evaluation (computing all 72 derivative quantities for 24 fields) accounts for roughly 40% of computation time. Replacing explicit 7-point 6th-order stencils with compact FD stencils reduces the per-derivative flop count (fewer explicit multiplications) at the cost of one tridiagonal solve per direction per field. The net benefit depends on the ratio of solve cost to stencil evaluation cost, which favors CFD at high accuracy requirements (where explicit schemes need wider stencils).

### 7.2 Higher-Dimensional Spectral Analysis

One open topic from E1: we are interested in extending our spectral analysis to higher dimensions, analogous to Lele 1992 Figure 4 and Chu & Fan 1998 Figure 1, which show **polar plots of phase speed anisotropy** for various schemes.

These plots are constructed as follows:

1. For a 2D grid with spacing $h$, consider a wave propagating at angle $\theta$ to the $x$-axis with wavenumber magnitude $\omega$.
2. The phase speed in direction $\theta$ is $c_\theta = \omega / k_\theta$ where $k_\theta$ is the numerical wavenumber in direction $\theta$.
3. In a Cartesian-stencil scheme, the numerical wavenumber at angle $\theta$ is $k_\theta(\omega) = \sqrt{\bar{\kappa}^2(\omega\cos\theta) + \bar{\kappa}^2(\omega\sin\theta)}$ where $\bar{\kappa}$ is the 1D modified wavenumber.
4. Plot $|c_\theta / c_{\text{exact}}|$ vs angle $\theta$ in polar coordinates; a perfectly isotropic scheme gives a circle.

For CCD, Chu & Fan Figure 1c shows nearly circular polar plots, indicating excellent isotropy. This is directly relevant to our operators since anisotropic phase errors in BSSN can cause non-physical features in gravitational waveforms.

Implementing this analysis in Mathematica is straightforward once we have the 1D modified wavenumber $\bar{\kappa}(\kappa)$ as a function. It should be added as a diagnostic section in `CFD Eigenvalue Stability.nb` or a new notebook.

---

## 8. Open Questions

1. **Convergence order at reduced $r$.** Why do A6 2nd-deriv operators stable at $r \leq 0.70\pi$ converge at ~2nd order in evolutions? Diagnose interior vs. boundary contribution (Track A, §A4).

2. **B-scheme failure mechanism.** Is the low convergence order of B-schemes an interior accuracy bug, a boundary closure inconsistency, or a ghost zone incompatibility? (Track B).

3. **CCD boundary polynomial.** How should the Kim24 polynomial framework be adapted when both $f'_j$ and $f''_j$ conditions are imposed? Does the degree-13 polynomial suffice, or must it be higher degree? (Track C, step 2).

4. **Outlier eigenvalues and the non-Dirichlet probe.** Does $r = 0.7$ in the eigenvalue test change the classification of operators with outlier eigenvalues? (Track A, §A2).

5. **Zhou-Zhang weight: minima shift.** Does using $\xi(\kappa) = D(\kappa)^2$ as the integration weight shift the locations of optimal interior parameters significantly relative to the unweighted case? (Track A, §A5).

6. **CCD memory feasibility.** At 24 BSSN fields × 2 (for $f'$ and $f''$) × 3 spatial directions = 144 derivative arrays per grid point, is the memory overhead compatible with Dendro-GR's GPU implementation? Verify before committing to CCD.

7. **S6 derivation typo.** There is a suspected typo in the S6 (heptadiagonal LHS) derivation. Diagnosing this is low priority and deferred.

8. **8th-order formal accuracy.** Deferred per v1.1 §4 B2 response. Wider stencil + more free parameters at 6th order is preferred over higher formal order at the same stencil width.

---

## 9. Step-by-Step Next Steps Summary

For each track, the immediately actionable work is:

**Track A (start now):**
1. Fix `yHat1List[[5]]` index bug.
2. Extend $\hat{y}$ scan range to $(-6, 6)$ in boundary notebook.
3. Run full operator enumeration at $r = 0.70\pi$ for A6 2nd-deriv (not just 3 operators); test all stable operators for convergence.
4. Diagnose convergence order: periodic BCs first to isolate interior behavior.
5. Implement non-Dirichlet probe ($r = 0.7$) and compare with current results.
6. Once convergence anomaly is understood, run hybrid Kim/Tyler scan for A6 2nd-deriv interior.
7. Prototype Deshpande KKT on A6 2nd-deriv interior as a closed-form alternative to `NMinimize`.

**Track B (start in parallel with A):**
1. Compare B6 interior modified wavenumber against A6 for 2nd derivatives.
2. Verify Taylor constraint derivation for 9-point RHS B-schemes.
3. Run periodic-BC convergence test for B6 2nd-deriv interior scheme.

**Track C (begin after A gives clear picture):**
1. Read Hixon 2017 closely for $\sigma_1$, $\sigma_2$ parameterization.
2. Derive joint CCD interior modified wavenumber in Mathematica.
3. Formulate CCD boundary polynomial constraint structure.
4. Implement new CCD notebooks once derivations are in hand.
5. Assess memory feasibility for Dendro-GR integration.

---

## 10. Additional Literature (from Web Search)

The following papers were identified in a supplementary literature search as potentially relevant. They are provided for reference; detailed reading is not required at this stage.

**Directly relevant:**

- **Daszuta et al.** (2024). "Spectrally-Tuned Compact Finite-Difference Schemes with Domain Decomposition and Applications to Numerical Relativity." *Journal of Computational Physics*. DOI: 10.1016/j.jcp.2024.002079.
  The **first published application of compact FD to numerical relativity**, using the Z4c system in GR-Athena++. Introduces a domain-decomposition strategy with approximate dispersion-relation preservation. Reports at least an order-of-magnitude reduction in phase error compared to standard explicit FD. **This paper should be read before beginning Track C**, as it directly addresses many of the implementation questions we face.

- **Mattsson & Ljungberg Rydin** (2022). "Implicit Summation by Parts Operators for Finite Difference Approximations of First and Second Derivatives." *Journal of Computational Physics*, 111743.
  Develops implicit (compact) SBP (summation-by-parts) operators with penalty boundary conditions (SAT method), up to 8th-order convergence. The SBP-SAT framework provides rigorous stability proofs for initial-boundary-value problems. A different approach from Kim24 but potentially complementary for understanding boundary stability.

- **Wang et al.** (2018). "A New Central Compact Finite Difference Scheme with High Spectral Resolution for Acoustic Wave Equation." *Journal of Computational Physics*, 366, 191–206.
  Combines cell-node and cell-center values to approximate second-order spatial derivatives with spectral-like resolution. Related to CCD but distinct framework.

- **Pan et al.** (2025). "Combined Compact Finite Difference Schemes for 2-D Acoustic Wave Propagation." *Journal of Applied Mathematics and Computing*, 71, 3009–3035.
  8th-order-in-space CCD for 2D acoustics. Relevant context for CCD accuracy in higher dimensions.

**For background/context:**

- **Brady & Livescu** (2019). *Computers & Fluids*, 183, 84–101. *(Now in Papers directory.)* High-order stable and conservative boundary schemes for explicit and compact FD at orders 4, 6, 8. See §6.9 of this document for discussion.

- **Pradhan et al.** (2015). "Construction, Analysis and Application of Coupled Compact Difference Scheme in Computational Acoustics and Fluid Flow Problems." *Communications in Computational Physics*, 18(4). Simultaneous computation of 1st, 2nd, and 4th derivatives via CCD; DRP performance analysis.

---

## 11. Key References (Updated Paths)

All papers are now in `C:\Users\natha\Documents\School\Research\Hirschmann\ClaudeCFD\Papers\`.

- `Kim07.pdf` — Kim (2007), interior optimization via spectral error functional
- `Kim24.pdf` — Kim (2024), stable boundary closure via extrapolating polynomial
- `Lele 1992.pdf` — Lele (1992), foundational compact finite differences; Table V and Figure 4
- `Jonathan Tyler Thesis.pdf` — Tyler BYU Master's thesis; total-error spectral optimization
- `Chen 2023.pdf` — Chen et al. (2023), integration-band cautionary tale for 2nd derivatives
- `Hixon CCD 2017.pdf` — Hixon (2017), CCD with $\sigma_1$, $\sigma_2$ optimization and b5 boundary
- `Chu Fan CCD 1998.pdf` — Chu & Fan (1998), foundational CCD; joint $f'$/$f''$ unknowns
- `Chu Fan 2000.pdf` — Chu & Fan (2000), staggered CCD (SCCD)
- `Deshpande 2019.pdf` — Deshpande et al. (2021, preprint 2019), KKT framework for optimized compact FD
- `Zhou Zhang 2011.pdf` — Zhou & Zhang (2011), prefactored optimized compact FD for 2nd derivatives; rational weight $\xi(\kappa)$
- `Brady Livescu 2019.pdf` — Brady & Livescu (2019), stable conservative boundary schemes
- `Dendro-GR 2022.pdf` — Dendro-GR GPU paper (2022)
- `Dendro-GR 2023.pdf` — Dendro-GR binary black hole paper (2023)
- `Gurarslan CCD 2021.pdf` — Gurarslan (2021), CCD for advection-diffusion with variable parameters
- `Zhang 2016.pdf` — Zhang (2016), splitting algorithms for high-bandwidth FD solves

**Mathematica notebooks:** `C:\Users\natha\Documents\School\Research\Hirschmann\HirschmannGroupBYU\Notes\MathematicaNotebooks\`

- `CFD Interior Optimization.nb`
- `CFD Boundary Optimization.nb`
- `CFD Eigenvalue Stability.nb`

**Additional literature to acquire:** Daszuta et al. 2024 (*J. Comput. Phys.*) should be obtained and added to the Papers directory before beginning Track C.
