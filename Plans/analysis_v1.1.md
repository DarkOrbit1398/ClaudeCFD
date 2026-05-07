# Analysis: Path to 6th-order Compact Finite Difference 2nd-Derivative Operators

Date: 2026-04-30
Source: synthesis of Kim07, Kim24, broader CFD literature, and the current state of the three Mathematica notebooks (`CFD Interior Optimization.nb`, `CFD Boundary Optimization.nb`, `CFD Eigenvalue Stability.nb`).

This document is a **proposal for review**, not a decided plan.

---

## 1. Background recap

### 1.1 Kim07 (interior optimization)

- 4th-order central pentadiagonal-LHS / 7-point-RHS first-derivative scheme.
- 5 unknowns $\{\alpha, \beta, a_1, a_2, a_3\}$; 2 consumed by formal Taylor constraints (Eqs. 2–3 of Kim07); **3 free parameters** for spectral optimization.
- Spectral-error functional (Kim07 Eq. 7):

  $$E = \int_0^{r}\!\Big\{(1+\delta)\kappa\big[1+2\alpha\cos\kappa + 2\beta\cos 2\kappa\big] - 2\sum_{m=1}^{3} a_m\sin(m\kappa)\Big\}^{2}\!\!\Big(\frac{\kappa}{r}\Big)^{n}\,d\kappa$$

  with $r \approx 0.85\pi$, $\delta = -2.33\times 10^{-4}$, $n = 10$.
- Boundary closure in Kim07 used a polynomial+trigonometric spline with control variables tuned by Newton-Raphson. This is the part Kim24 supersedes.
- **First derivatives only.** Kim07 says nothing about $f''$.

### 1.2 Kim24 (boundary closure)

- Single degree-13 polynomial $P(\hat{x}) = \sum_{m=0}^{13} p_m \hat{x}^m$, determined by 14 constraints:
  - 7 function-value matches: $P(j) = f_j$ for $j = 0, \ldots, 6$
  - 6 derivative matches: $P'(j) = f'_j$ for $j = 0, \ldots, 5$
  - 1 free constraint: $P^{(n)}(\hat{y}) = 0$ with $n = 8$
- $\hat{y}$ is a **single shared free parameter** across all 3 boundary rows in the original Kim24 paper.

> NATE COMMENT: Where does Kim indicate that he used the same value of $\hat{y}$ for all rows in his paper? I don't recall him ever specifying what was done with $\hat{y}$, but at the end of Section 2.4 he states "once the free variables $\hat{y}$ are determined by the optimisation." The plural *variables* indicates to me that he may be using a different $\hat{y}$ for each boundary row just as we are.

- Three boundary rows ($i = 0, 1, 2$) at each end; LHS retains pentadiagonal structure, RHS uses 7 points.
- **Optimization objective is different from Kim07**: physical-space, multi-grid $L_\infty$ error on $\exp(\sin x)$ across 8 grid resolutions, not spectral error.

> NATE COMMENT: Sure, he is using $L_\infty$ error optimization rather than spectral error optimization for the boundary nodes. However, at the beginning of Section 2.2, he states "while the fourth-order interior scheme optimised by Kim [9] for maximised resolution performance is directly adopted." The reference Kim [9] in the Kim24 paper is a direct reference to the Kim07 paper, indicating that he is still using spectral error optimization for the interior nodes. What differences, benefits, or drawbacks might we see from both of these different methods of error optimization? For me personally I have thought the spectral error minimization is better because the spectrum is specific to an operator, while an $L_\infty$ error minimization is specific to a test function. Thus if we want to apply the same CFD operator to different test functions (testing in easier systems like Maxwell's equations and ultimately using the same operator for the BSSN equations), we would prefer spectral error minimization.

- Stability constraints enforced **during optimization**:
  - Non-Dirichlet weighted boundary probe ($r = 0.7$ in Kim24 Eq. 17) — exposes instabilities Dirichlet analysis hides.

  > NATE COMMENT: I don't really understand what a non-Dirichlet weighted boundary probe is, and I haven't had the time yet to understand what Kim is doing with this $r = 0.7$ value in his paper. I am hoping to come back to this, but if I don't, please explain better what this is and what its purpose is in the context of the comapct finite difference problem or in the context of numerical methods more generally.

  - $\text{Re}(\lambda_i) \geq 0$ (eigenvalue criterion).

  > NATE COMMENT: In our sign convention, our stability criterion is that $\text{Re}(\lambda_i) \leq 0$, which is opposite Kim's sign convention. I think it is also important to note that this is the stability criterion for first derivatives. For second derivatives, a different stability criterion must hold.

  - $|R(-\text{CFL}\cdot\lambda_i)| < 1$ at CFL = 1.047 (RK4 polynomial bound).
- Solver: `scipy.optimize.differential_evolution`.
- **First derivatives only.** Kim24 says nothing about $f''$.

### 1.3 Our extensions to date

- Mathematica notebooks reproduce Kim07 interior and Kim24 boundary frameworks for first derivatives across A4/B4/C4/A6/B6/C6/A8 schemes.
- Extended to second derivatives at A4 and A6 only.

> NATE COMMENT: This is not true, I have extended the second derivatives for B4/C4/B6/C6 as well, but progress in the notebook was lost due to Git overriding my file with a previous version. I will fix this at some point, hopefully soon. I also need to check if it's this way in the interior and boundary notebooks or just the interior notebook (I thought it was just that one that got overridden).

- Deviated from Kim24 by using **one $\hat{y}$ per boundary row** (yHat0, yHat1, yHat2) rather than a single shared $\hat{y}$ — confirmed working for first derivatives.

> NATE COMMENT: Again, is this a deviation from Kim's method, or was his paper misinterpreted?

---

## 2. Diagnosis: why current 2nd-derivative operators plateau at "explicit-FD-equivalent" accuracy

Smoking-gun number from `CFD Interior Optimization.nb` (stored output, A6 second-deriv):

| metric | Nate (our optimized) | Tyler/Kim baseline |
|---|---|---|
| low-wavenumber error | **4.5e-7** | 4.7e-6 |
| total error | 1.6e-3 | **2.65e-4** |

Our optimization wins low-$\kappa$ by 10× but loses total error by 6×. Classic over-fitting of the low-wavenumber regime at the cost of high-$\kappa$ behavior. Five contributing causes follow.

> NATE COMMENT: I believe it got deleted with the version overriding mentioned earlier, but there should be similar tests for all the other operators, not just A6 second derivatives. We have understood for some time that Kim's approach minimizes the low-wavenumber spectral error, while Tyler's method minimizes the total spectral error. This is why I developed hybrid interior schemes that take some free parameters to be "Kim-like" and some to be "Tyler-like" in order to meet in some middle ground between these two error windows. Again, the hybrid schemes may have been lost to the version overriding, I will have to double check and maybe code them back up again.

### 2.1 Lost free parameter from the consistency constraint

Second derivatives require an extra Taylor constraint that first derivatives do not: the row sum of Q must vanish (encoded in the notebook as the `−2s` diagonal term, with $s = a_1 + a_2 + a_3$).

| scheme | LHS | RHS | free params |
|---|---|---|---|
| 1st-deriv A6 | pentadiagonal | 7-point | **3** ($a_1, a_2, a_3$) |
| 2nd-deriv A6 | pentadiagonal | 7-point | **2** ($a_2, a_3$) |

One fewer optimization knob at the same stencil width. Cannot match what 1st-deriv achieves.

> NATE COMMENT: This is not true. The same constraint arises in the first derivatives that the right hand side (RHS) parameters sum to zero (since we have $-a_i$ on the left of the main diagonal and $+a_i$ on the right of the main diagonal, so it is purely antisymmetric). Of course that is for the interior, but one finds that the sum of the parameters in each boundary row also sums to zero. I think this is the zeroth order Taylor constraint(?) but I forget exactly.

### 2.2 Error functional weight is wrong for 2nd derivatives

Kim07's $(\kappa/r)^n$ weight with $n = 10$ is appropriate for $f'$ where the modified-wavenumber target is **linear** in $\kappa$. For $f''$ the target is $\kappa^2$, and the integrand $((1+\delta)\kappa^2 D - N)^2$ is already $\kappa^4$-scaled before the weight is applied. Net effect: the high-$\kappa$ tail dominates the cost, low-$\kappa$ is under-served, and the integral runs to $r \approx 0.85\pi$ where Chen et al. 2023 demonstrate that 2nd-derivative optimization develops harmful overshoot in the modified wavenumber.

> NATE COMMENT: We had discarded the parameters $n$ and $\delta$ from the Kim07 error optimization, thinking they were just additional "bells and whistles" that seemed only to make a negligible difference to the resulting interior coefficients when we tested it, and it seemed to complicate the spectral error function more than was necessary. From this comment it sounds like they do play a role, but I do not understand what their role is? I will read the Chen 2023 paper to see if it clears things up but I would like an explanation of why Kim includes $n$ and $\delta$ in his spectral error function, and why that would matter for first vs. second derivatives.

Chen et al. (2023) recommend integrating only to $r \approx 0.5\pi$–$0.75\pi$ for 2nd-derivative spectral optimization.

> NATE COMMENT: We chose $r \sim 0.85 \pi$ rather arbitrarily as well; it just seemed like a choice that worked well. Is there a more systematic approach to choosing a good value? Is it problem-dependent based on how high of wavenumbers/frequencies we would like to resolve well? Perhaps the Chen 2023 paper addresses this as well.

### 2.3 Boundary $\hat{y}$ is hand-picked, not optimized

`CFD Boundary Optimization.nb` runs a 1D grid scan `Range[0, 5, 0.03]` per node, plots, and the user manually selects local minima. With three boundary rows (per-row $\hat{y}$), the search space is 3D and only ~60 combinations are tested combinatorially in the Eigenvalue notebook.

> NATE COMMENT: The number of operators tested depends on the number of minima found for $\hat{y}$ in each boundary node; some only have about 60 combinations, but others have hundreds.

Implications:
- Optima between grid points are missed.

> NATE COMMENT: I chose to only test a discrete set of values for each $\hat{y}$ since the spectral error function involved a rather complicated integral that had to be evaluated numerically yet still took a long time. I played around with different resolutions and ranges, and going from $0$ to $5$ in steps of $0.03$ seemed to produce a smooth error function with clear minima; for some operators I found it helpful to extend the range into the negatives, but not very far (for example, there were some minima with values such as $\hat{y} \sim -0.1$). I don't believe that we are missing minima with our current procedure.

- Stability (eigenvalue location) is not enforced **during** boundary optimization — it's verified after the fact by trimming. Accuracy and stability fight each other rather than co-optimizing.

> NATE COMMENT: This is true. We are first optimizing for accuracy and then checking for stability after. As far as I have seen in the literature, there isn't really a good way to develop operators for stability rather than just trying a bunch and checking which ones are stable (such as Brady and Livescu's approach where they ran a machine learning algorithm to find coefficients that produced stable operators). It may be good to have you read Brady and Livescu's paper to be familiar with their approach. If there are ways to optimize accuracy and stability at the same time, I would be curious to know what they are, but that seems implausible.

### 2.4 Boundary objective doesn't match Kim24's update

Kim24 abandoned spectral error as the boundary objective in favor of a **physical-space, multi-grid $L_\infty$ error**. Our boundary notebook still optimizes the spectral function. For first derivatives this works (we have stable, accurate operators). For second derivatives the spectral optimum may be more fragile.

> NATE COMMENT: What would cause the spectral error function for second derivatives to be more fragile than first derivatives? I know that even finite difference schemes have a different parity/symmetry structure than odd finite difference schemes so maybe that comes into play with the spectral error function.

### 2.5 Stability probe may still be Dirichlet

Kim24's whole point is that Dirichlet eigenvalue analysis masks instabilities the operator actually has when used in real problems. We need to confirm `CFD Eigenvalue Stability.nb` uses the $r = 0.7$ non-Dirichlet probe — if it's Dirichlet-only, we may be calling unstable operators "stable."

> NATE COMMENT: Again, I need to review Kim's paper to understand what this non-Dirichlet probe means and how it could mask instabilities. There have only been a few operators which seemed to be stable according to their eigenvalue spectra but were unstable when tested in evolutions like Maxwell's equations. These have all had features in common in the eigenvalue spectra which we have called "outlier eigenvalues," and they still have nonpositive real part but they have larger than normal imaginary parts, so they still sit in the left side of the complex plane, but they deviate from the normal shape/pattern that an eigenvalue spectra creates.

---

## 3. Notebook state summary

| notebook | 1st deriv | 2nd deriv |
|---|---|---|
| Interior Optimization | T2, A4, B4, C4, A6, B6, C6, A8 + spectral hybrids | A4, A6 only |
| Boundary Optimization | A4, B4, C4, A6, B6, C6, S6 (stuck), A8 | A4, B4, C4, A6, B6, C6 (no S6, no A8) |
| Eigenvalue Stability | A4, B4, C4, A6, B6, C6, A8 | A4, B4, C4, A6, B6, C6 (no A8) |

Known issues flagged by the notebook inventory:

- 2nd-deriv P6(A6) `Output Coefficients` cell has a live `Part::partw` error: `yHat1List[[5]]` indexes a length-4 list. Needs fixing before re-running.

> NATE COMMENT: I am pretty sure this is just an issue where I was working on the code, copied some values from another section as a placeholder, and then forgot to change it in another location. I will fix it.

- "Build Modified Matrices (needs cleaning badly)" labels in multiple boundary subsections.

> NATE COMMENT: I just think the code in that block is a bit messy and hard to read. It functions correctly and can be explained, but I just wanted to optimize it a little bit to be more readable. This is a low priority issue, but any tips on how to make the code in this section more readable would be appreciated.

- S6 1st-deriv flagged "ALMOST DONE, JUST STUCK ON THE EIGENVALUES."

> NATE COMMENT: The A/B/C operators that we have developed are all pentadiagonal on the left hand side (LHS). I was developing a septadiagonal LHS 6th order scheme for first derivatives, I had made a lot of progress but the eigenvalue spectra weren't what I anticipated. I think there might have been an issue in my derivation and so I paused and never got around to finishing it. I could fix this sometime if we decide that pursuing septadiagonal (or heptadiagonal, just 7-point on the LHS) schemes is worthwhile.

---

## 4. Proposed plan: three tracks

### Track A — Tighten what we already have (1–2 weeks, high-confidence gains)

Low-risk, high-leverage edits to the existing pentadiagonal-LHS / 7-point-RHS A6 second-derivative scheme.

**A1. Narrow the spectral-error integration band for 2nd derivs.** Drop $r$ from $\sim 0.85\pi$ to $0.65\pi$–$0.75\pi$ for the second-derivative error functional only (Chen et al. 2023). One-line change. Re-run A6 2nd-deriv optimization and compare total-error metric.

> NATE COMMENT: This could be worth doing. I think we ought to read the Chen paper to see what benefits they saw from a lower cutoff wavenumber $r$ or we ought to address the question above of how to systematically choose a good value of $r$, but I agree that this is a rather easy fix and would be quick to test.

**A2. Reweight the error functional for $\kappa^2$ targets.** Replace $(\kappa/r)^{10}$ with a weight that matches the $f''$ target. Two candidates:

- Zhou-Zhang weight: $\xi(\kappa) = (1 + 2\alpha\cos\kappa + 2\beta\cos 2\kappa)^2$ — makes the integrand rational, cleaner objective.
- Relative-error weight: $\big((\kappa^2 D - N)/\kappa^2\big)^2$ — penalizes fractional error, which is what we actually care about.

> NATE COMMENT: These seem like things which could be worth exploring, but I don't understand their purpose or motivation. I understand at least why making the integrand rational would be helpful as that's the reason we used $(\kappa D - N)$ instead of $(\kappa - \tfrac{N}{D})$, but with the Zhou-Zhang weight, how do we know that this weight wouldn't change the behavior of the spectral error function, or more specifically where the minima are?

**A3. Fix the `yHat1List[[5]]` index bug** in 2nd-deriv P6(A6) before any re-run.

> NATE COMMENT: I'll fix this. Again, it's not really an issue.

**A4. Verify the eigenvalue probe.** Confirm `CFD Eigenvalue Stability.nb` uses the non-Dirichlet $r = 0.7$ probe (Kim24 Eq. 17). If not, switch.

> NATE COMMENT: We are *not* using the non-Dirichlet probe. I would like to understand its utility more before we consider implementing it as it seems more problem-specific and less operator general.

**A5. Replace manual yHat scan with `NMinimize`** over the joint (yHat0, yHat1, yHat2) space subject to `Max[Re[Eigenvalues[M]]] ≤ 0` as an inequality constraint. Mathematica supports this directly.

> NATE COMMENT: I will try to see if I can implement `NMinimize` into the boundary spectral error optimization. I did not know about that command previously, so I will check it out.

### Track B — Widen the stencil (2–4 weeks)

To recover the missing free parameter relative to 1st-deriv A6.

**B1. Heptadiagonal RHS (`nR = 4`, 9-point RHS) at pentadiagonal LHS.** Adds $a_4$ to the RHS, giving $\{a_2, a_3, a_4\}$ free with $\{a_1\}$ Taylor-constrained. Cost: one more MPI ghost-zone layer. Most natural extension and matches how 1st-deriv A8 widens.

> NATE COMMENT: We have tried this approach already. These are the so-called "B" schemes in the Mathematica notebooks (B4/B6). However, we saw very poor performance of these operators in testing. It is possible that there are issues in the code in comparison with the A/C schemes, but I am unsure and that's something we'd still like to understand at some point.

**B2. On 8th order specifically.** With $(n_L = 2, n_R = 4)$ and 8th-order Taylor target, you consume 4 of 5 RHS coefficients, leaving only 1 free — not better than A6. So 8th order **as a target order** isn't the lever. The right framing: **wider stencil + same 6th-order target + more free parameters for spectral optimization**. This is exactly Kim's philosophy applied at a wider stencil.

Lele 1992 Table V already shows the resolution-efficiency gain from 6th to 8th *formal* order is small. Don't pursue 8th formal order. Pursue wider stencils with more free parameters at 6th formal order.

> NATE COMMENT: We will hold off on 8th formal order for now. If we pursue 8th formal order, it seems like it will be necessary to move to a heptadiagonal LHS rather than a pentadiagonal LHS. There was an attempt of this with the S6 operator in the Mathematica notebook, but it is unfinished and can be left for another time.

### Track C — Pivot to Combined Compact Difference (4–8 weeks, biggest payoff)

The literature's strongest answer to the 2nd-derivative accuracy plateau.

**Why CCD.** Chu & Fan (1998) treat $f'$ and $f''$ as joint unknowns at every node, coupling them via two compact relations derived from a single degree-6 Hermite polynomial. Leading-error coefficients (per-step $h^p$):

| scheme | $f'$ leading error | $f''$ leading error |
|---|---|---|
| 6th explicit central | ~36/7! | ~72/8! |
| Lele 6th tridiagonal | ~4/7! | ~17/8! |
| Lele 6th pentadiagonal | ~5/7! | ~28/8! |
| **CCD** | **~0.87/7!** | **~2/8!** |

CCD's $f''$ is **~8× more accurate** than the 6th-order tridiagonal compact, on a 3-point stencil. Cost: ~2× a single tridiagonal solve (twin-tridiagonal block elimination), 2× memory per field.

**Directly applicable template.** Hixon 2017 (`HixonCCDOptPaper.pdf`) takes Chu-Fan CCD and adds:

- Two free parameters $\sigma_1, \sigma_2$ in the leading $f'$ and $f''$ error terms, optimized for spectral resolution (Kim07's idea applied to CCD).
- Boundary closures b1/b4/**b5**, where b5 matches the leading-error term to recover full design order.
- Eigenvalue-stability verification of the resulting joint operator.

**What porting looks like.** Two new notebooks paralleling the existing structure:

- `CCD Interior Optimization.nb` — Kim07-style spectral-error integral on the $f''$ modified wavenumber subject to the joint $f'$-$f''$ coupling.
- `CCD Boundary Optimization.nb` — extend Kim24's polynomial-extrapolation trick to the **coupled** (f', f'') unknown vector. The single polynomial would need to supply consistent $f, f', f''$ extrapolations. This is genuinely new derivation work, but the Kim24 framework transfers cleanly: degree-13 polynomial, $n = 8$ zero-derivative constraint, free $\hat{y}$ per boundary row, joint stability optimization.

**Verdict.** If we want second-derivative accuracy that meaningfully exceeds explicit FD, CCD is the most promising direction in the literature. Tracks A+B may close half the gap to first-derivative-quality results within the existing scheme framework; CCD likely closes the rest.

> NATE COMMENT: We are not familiar with the CCD approach, and I would like to read the Chu Fan and Hixon papers to become familiar with this approach before we move forward and try it. A more clear summary of the benefits and drawbacks of the CCD approach would be helpful information.

---

## 5. Other ideas considered and not recommended

| Idea | Source | Verdict |
|---|---|---|
| Apply $D^1$ twice to get $D^2$ | Sari & Gürarslan | Worse than dedicated $D^2$. Skip. |
| Spotz-Carey HOC | Spotz-Carey 1996 | Substitutes the PDE into the truncation error. Not a generic $D^2$ operator. Skip. |
| Ciment-Leventhal $(I+\delta^2/12)^{-1}\delta^2$ | Ciment-Leventhal 1975 | Wave-equation-specific. Skip. |
| Lee-Liu-Sun CCD2 (mixed derivs) | 2014 | Premature for our 1D needs. |
| Central Compact Schemes (Liu/Shu) and UTEP CHVM | 2013, 2022 | Double per-point DOF; painful with ~25 BSSN fields. Skip unless 1st-deriv resolution becomes a separate bottleneck. |
| Splitting algorithms for high-bandwidth solves | Zhang 2016 | Cost reduction, not accuracy gain. Useful only if we exceed pentadiagonal. |
| Zhou-Zhang prefactorization | Zhou-Zhang 2011 | Cost / MPI-communication win, not accuracy. Worth knowing for production. |
| Deshpande KKT framework for analytic optimization | Deshpande 2021 | Closed-form alternative to numerical optimization. Could replace `NMinimize` cleanly. **Worth considering as a tooling improvement.** |

> NATE COMMENT: We have evidence that applying a first derivative twice is worse than a second derivative once, so we are happy to skip that approach. I am unfamiliar with the rest of the techniques, so I cannot speak to those. However, if we are going to rule them out, please provide more information on why we ought to skip them; the brief sentences in the verdict here are not enough.

---

## 6. Recommended sequencing

1. **Start Track A immediately** — low cost, diagnostic. Tells us how much of the gap is over-fitting vs. structural.
2. **Start Track C derivation in parallel** — the boundary-polynomial extension to coupled (f', f'') is the hard part and worth several weeks of derivation work regardless of Track A's results.
3. **Decide Track B based on Track A's outcome** — if A doesn't close the gap, widen the stencil before going to CCD.
4. **Defer 8th-order formal target.** Use any wider stencil for free parameters, not for additional Taylor constraints.

> NATE COMMENT: Before we proceed with any specific next steps, please address all the comments provided and study the additional material provided so that we have a better-informed second draft of this proposal.

---

## 7. Open questions

These need answers before any code is written.

1. **Eigenvalue probe**: does `CFD Eigenvalue Stability.nb` already use Kim24's non-Dirichlet $r = 0.7$ probe, or only Dirichlet? The notebook agent didn't surface this explicitly.

> NATE COMMENT: We are *not* using Kim24's non-Dirichlet probe. We want to optimize our CFD operators more generally, not to a specific test function as in Kim24's approach with the $L_\infty$ error optimization. This is why we are opting to extend and generalize Kim07's spectral error optimization approach to the boundaries, incorporating the extrapolation polynomial technique from Kim24.

2. **Per-row vs shared $\hat{y}$**: project.md says we have one $\hat{y}$ per boundary row (extension of Kim24's shared single $\hat{y}$). Confirmed working for first derivatives, or still experimental?

> NATE COMMENT: We have confirmed that using a different value of $\hat{y}$ for each boundary node is okay. Kim24's approach does not specify explicitly whether a single value of $\hat{y}$ is used or whether a different $\hat{y}$ value is used for each boundary node.

3. **MPI patch boundaries**: how are interior/exterior patch boundaries handled — same boundary closure as physical boundaries, or overlapping ghost zones with interior stencil? Affects how aggressively boundary closures matter.

> NATE COMMENT: For more information on how we handle boundaries, ghost zones, etc., please refer to the Dendro-GR papers which will be provided.

4. **Test functions**: Kim24 uses $\exp(\sin x)$. Numerical relativity has its own test cases (linearized waves on a black-hole background, etc.). Should the optimization objective use a relativity-relevant test function?

> NATE COMMENT: As mentioned above, we do *not* want to optimize the CFD operators to a specific test function, even if it is relevant to Numerical Relativity. We will develop the operators more generally and then test them in several different systems, including Maxwell's equations, the Nonlinear Sigma Model, Teukolsky waves, etc., but this should not yet be relevant to the development of the operators.

5. **Computational budget**: are we willing to take a 2× memory hit per field for CCD, given ~25 evolved BSSN/Z4c fields per point?

> NATE COMMENT: Probably not, but we need to understand the benefits and drawbacks of the CCD method more clearly before making a definite conclusion.

6. **Tyler/Tyler-BYU**: the notebook references "Tyler" coefficients and the literature directory has `Tyler-BYU/`. Prior BYU work we should also be reading?

> NATE COMMENT: Yes, please review Jonathan Tyler's BYU Master's Thesis, which we commonly refer to as the Tyler paper.

> NATE COMMENT: Please review all comments made in this document. Please review the additional literature provided, including the Brady and Livescu paper, the Tyler paper, and the Dendro-GR papers. Please make the necessary adjustments to this document to address all comments and new information. Please ask us questions if anything remains unclear or if more information is needed.

---

## 8. Key references (paths)

- `/Users/dave/Research/compact_finite_difference/Kim/Kim-optimization.pdf` — Kim07
- `/Users/dave/Research/compact_finite_difference/Kim/kim-stable-boundary-closure.pdf` — Kim24
- `/Users/dave/Research/compact_finite_difference/CCD/threepoint_ccd_scheme.pdf` — Chu & Fan 1998 (foundational CCD)
- `/Users/dave/Research/compact_finite_difference/HixonCCDOptPaper.pdf` — Hixon 2017 (CCD optimization + b5 boundary, the directly-applicable template)
- `/Users/dave/Research/compact_finite_difference/Zhou-CFD-Wave-Eq-2nd-Deriv.pdf` — prefactorization for 2nd derivs
- `/Users/dave/Research/compact_finite_difference/Deshpande-generate-optimized-central-finite-difference.pdf` — closed-form KKT optimization with stability
- `/Users/dave/Research/compact_finite_difference/2nd-derivs-remotesensing-15-00604-v2.pdf` — Chen 2023, integration-band cautionary tale
- `/Users/dave/Research/compact_finite_difference/Lele_1992JCP.pdf` — Lele 1992, foundational
