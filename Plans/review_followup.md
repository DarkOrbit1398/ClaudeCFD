# Follow-up to analysis_v1.1.md review

Date: 2026-05-01
Purpose: capture all open items from the v1.1 NATE COMMENT review so that nothing is lost between this conversation and the v1.2 redraft.

Each item is numbered. Please respond inline under each (e.g. with a "RESPONSE:" block) — no need to keep the numbering tidy.

---

## A. Corrections to v1.1 that I plan to make in v1.2

These are points where your comments showed v1.1 was wrong or incomplete. Confirm the correction is right, or push back if I have it wrong.

### A1. Free-parameter accounting (v1.1 §2.1) was wrong

Your antisymmetry argument is correct: 1st-deriv RHS coefficients also satisfy a row-sum-zero constraint (zeroth-order consistency); it is just hidden by the antisymmetric form. The "smoking gun" of fewer free parameters in 2nd-deriv A6 needs to be re-derived from scratch against the actual notebook definitions. I will redo the parameter count carefully in v1.2.

> RESPONSE: The first part is correct. I don't know that any parameter counting leads to a "smoking gun" here, but please double check the parameter counting in v1.2.

### A2. Kim24 boundary objective — I overstated the spectral→L∞ "switch" (v1.1 §1.2)

Kim24 keeps Kim07 *spectral* error for the **interior** (the Kim [9] reference) and only uses physical-space $L_\infty$ for **boundary closure tuning**. Your point about test-function-independence shifts the v1.2 recommendation: continuing spectral-error optimization at the boundary (with Kim24's polynomial as the technique) is principled, not a mistake. v1.2 will reflect this.

> RESPONSE: This is correct, the choice to use the spectral error optimization in our approach is intentional.

### A3. Per-row vs shared $\hat{y}$ — calling this a "deviation" is unjustified (v1.1 §1.2, §1.3)

Kim24's Section 2.4 plural ("variables $\hat{y}$") is genuinely ambiguous. v1.2 will not assert that per-row $\hat{y}$ is a deviation; instead I will re-read the optimization tables in Kim24 and report what the paper's reported coefficients actually disambiguate.

> RESPONSE: Upon reviewing the Kim24 paper more closely, I think it is still rather ambiguous whether or not he uses a single value of $\hat{y}$ or a different value for each boundary node. You are welcome to try to disambiguate, but I think we should not get caught up on this detail, as we have already shown in our own work that using a separate $\hat{y}$ value for each boundary node works just fine.

### A4. Track B (heptadiagonal RHS) is not green field (v1.1 §4 B1)

B4/B6 already exist as 1st-deriv heptadiagonal-RHS schemes and have underperformed. v1.2 will reframe Track B as "diagnose why B-schemes lost," not "invent wider stencils." This depends on item B6 below.

> RESPONSE: This is correct. We are aware that wider stencils are an option which increase either our nominal order of accuracy or the number of free parameters available for optimizing the spectral resolution, but at increased computational cost due to larger grid sizes. It is correct that we have already investigated the B schemes, did not see the expected performance, and partially abandoned the idea. It would be helpful to investigate why this was the case.

### A5. The $(\kappa/r)^n$ weight and $\delta$ lift have been discarded in your notebooks (v1.1 §2.2)

My v1.1 diagnosis assumed Kim07's full functional was in use. v1.2 needs to re-derive the diagnosis from your actual functional. Please confirm scope (item B7 below).

> RESPONSE: This is correct. We are *not* using the full functional described in the Kim07 paper, but rather a similar version to it. I will provide below the specific functional we use, although I think this is more relevant to items B10 and B11 rather than item B7.

### A6. Hybrid (Kim + Tyler) interior schemes exist and were missing from v1.1

The "Nate vs Tyler" error gap was the *motivation* for hybrids, not a bug. v1.2 will discuss hybrids as a first-class option for $f''$ before any structural change.

> RESPONSE: This is true, it was never thought of as a bug to begin with. We understand that the two different methods for optimizing the spectral error will provide either better low-wavenumber error (in Kim's approach, where the spectral error is optimized up to some "cutoff wavenumber," which we have been generally using $r \sim 0.85 \pi$ but which we have discussed varying in other sections of this review) or better total-wavenumber error. (in Tyler's approach, which you will read. As a small note, Tyler discusses spectral schemes in Sections 3 and 4 of his thesis; Section 4 introduces his idea of including constraints of the form $\tilde{\omega}(\omega_i) = \omega_i$, which are the specific constraints we are referring to in this approach. In contrast, the constraints on the spectral function in Kim07's approach are seen in Eq. 8, where he sets $\tfrac{\partial E}{\partial \alpha} = 0$ and similar for the other free parameters. So then the hybrid approach allows for a combination of both types of these constraint equations.) Exploring the benefits/drawbacks of these differing spectral approaches is something we would like to do in v1.2.

---

## B. Information / files / clarifications I need before drafting v1.2

### B1. Path to Jonathan Tyler's BYU Master's thesis on disk — PROVIDED

Path: `/Users/dave/Research/compact_finite_difference/Tyler-BYU/BYU-Analysis and Implementation of High-Order Compact Finite Differen.pdf`

I will read this before drafting v1.2.

### B2. Path to the Brady & Livescu paper — PROVIDED

Path: `/Users/dave/Research/compact_finite_difference/Brady-Livescu/Brady_Livescu_CF19.pdf`

I will read this before drafting v1.2.

### B3. Dendro-GR papers — PROVIDED

Directory: `/Users/dave/Research/compact_finite_difference/Dendro-GR/`
Contents:
- `PhysRevD.111.124001.pdf`
- `PhysRevD.107.064035.pdf`
- `ics19.pdf`
- `IEEE Xplore Full-Text PDF.html`

Is there a specific one of these I should prioritize for the MPI patch / ghost-zone question, or read all of them?

> RESPONSE: I think the `PhysRevD.107.064035.pdf` paper would be a good primary reference, and the `IEEE Xplore Full-Text PDF.html` is another good choice. I think the other two aren't necessary. What I would like you to focus on in these papers again is understanding the Octree grid structure that Dendro-GR uses, as our goal is to ultimately implement the CFD operators for use in Dendro-GR. For example, I believe the blocks we use are $13 \times 13 \times 13$, and this ought to higlight to you that we have a limit on how large our stencil sizes for CFD operators are as well as the importance of our boundary closures performing well. You may find other information useful in the context of what computational costs come from replacing explicit finite differences with CFDs in the Dendro-GR code. However, please remember that reading these papers is intended for *context*, not for you to critique the methods of our implementation into the Dendro-GR code. Our primary goal in this project is still developing accurate and stable CFD schemes.

### B4. Chen 2023 — confirm I have the right paper

I assume Chen et al. 2023 = `/Users/dave/Research/compact_finite_difference/2nd-derivs-remotesensing-15-00604-v2.pdf`. Confirm or correct.

> RESPONSE: Correct, although the file path has changed (see Section E1 below).

### B5. Restoration of the version-overwritten notebook content

When B4/C4/B6/C6 second-deriv work and the hybrid interior schemes are restored, point me at the updated notebook(s) so I can re-inventory. Until then, treat the v1.1 notebook coverage table as out of date.

> RESPONSE: I will do this once the notebooks are updated. I don't know if "out of date" is the right way to treat the current notebook; I believe the current notebooks you have been provided are up-to-date in their methods, although they may only include some of the families of schemes we have tested, not all.

### B6. What "poor performance" of B-schemes actually means

Is the failure (a) spectral resolution at high $\kappa$, (b) eigenvalue stability, (c) time-domain instability in evolutions (e.g. Maxwell, NSM), (d) something else? This determines whether Track B is a derivation-bug hunt or a structural dead end.

> RESPONSE: I will have to double check this information since it has been a while, but as I recall it, the B-schemes had great spectral resolution (they had an additional free parameter over the A-schemes for spectral optimization). In terms of eigenvalue stability, many of the B schemes I tested were supposedly stable, I think some actually performed that way in evolution tests with the wave equation but some did not. The biggest issue I believe was that the order of convergence for the B-schemes was lower than the A-schemes and C-schemes.

### B7. Stability criterion for second derivatives

You noted the 1st-deriv criterion ($\text{Re}\,\lambda \leq 0$ in your sign convention) is the wrong question for $f''$. What criterion do you actually use? (For an operator used in $u_{tt} = c^2 u_{xx}$, my expectation is that the relevant eigenvalues should be real-negative or close to the negative real axis, but I want your specific convention.)

> RESPONSE: This is a good question. Let us formulate our CFD matrix system as $\mathbf{P} f' = \mathbf{Q} f$, or $f' = \mathbf{P}^{-1} \mathbf{Q} f$. Let $\mu$ be the eigenvalues of $\mathbf{P}^{-1} \mathbf{Q}$ such that $\mathbf{P}^{-1} \mathbf{Q} \vec{v} = \mu \vec{v}$ for an eigenvector $\vec{v}$. In this convention, then as you have mentioned, our stability criterion for second derivatives is that $\Re(\mu) \leq 0$ and $\Im(\mu) = 0$. However, we have also been using the convention that if we define $\lambda := \sqrt{\mu}$ such that $\mathbf{P}^{-1} \mathbf{Q} \vec{v} = \lambda^2 \vec{v}$, then our stability criterion for second derivatives is that $\Re(\lambda) = 0$ and $\Im(\lambda) \geq 0$ (hopefully it is not difficult to see that these two criteria are equivalent, it has mostly just been a matter of convention and you might see the line `eigs = (Eigenvalues[mat])^(1/2);` in the Mathematica code, indicating the use of the $\lambda$ convention as specified above).

### B8. The "outlier eigenvalue" pattern

Have you correlated outlier eigenvalues with (a) $\hat{y}$ values at a specific boundary row, (b) scheme family (A vs B vs C), (c) boundary-row index, (d) something else? Any systematic pattern would help me think about whether the non-Dirichlet probe targets them.

> RESPONSE: We have not correlated outlier eigenvalues with any specific structure in the CFD operators. There may be patterns that are not obvious, but the most I have noticed is that similar choices of $\hat{y}$ lead to similar CFD coefficients and thus similar eigenvalue spectra shapes (regardless of whether there are outliers or not). The outlier eigenvalues have not been holding up our analysis, so I would consider them a lower priority issue to understand. However, if there are simple tests we could do to investigate and understand their behavior, then we could do that at a later time.

### B9. Hybrid interior scheme structure

Even briefly: which parameters were "Kim-like" (low-$\kappa$ minimizing) vs "Tyler-like" (total-error minimizing)? Which scheme letter(s) do they live under? I can derive it back from the notebook once restored, but a sentence now saves a round trip.

> RESPONSE: I believe I addressed this in some depth in item A6 above. In the notebook `CFD Interior Optimization.nb`, go to the Section `Spectral Schemes`, Subsection `Tyler Spectral Schemes`, there are several Subsubsections which you may find helpful. `P4 (A4)`, `P6 (A6)`, and `P6 (B6)` have code to generate interior schemes where all free parameters are used as in the Tyler method, and `Hybrid A42K1T` and `Hybrid A41K2T` contain hybrid schemes (here the naming is supposed to be suggestive; `A4` meaning an A4 CFD scheme, which has 3 free interior parameters, and then `2K1T` or `1K2T` indicate that 2 free parameters are Kim-like and 1 is Tyler-like, or vice versa).

### B10. Confirm scope of $(\kappa/r)^n$ and $\delta$ removal

Was the $(\kappa/r)^n$ weight and $\delta$ lift discarded in (a) interior optimization only, (b) boundary optimization only, (c) both? And in (i) all schemes, (ii) just some — if just some, which?

> RESPONSE: It was discarded in the interior and boundary procedures, and for all schemes.

### B11. The notebook's actual current spectral-error functional

Tied to B10: please paste or point me at the exact integrand currently being minimized in `CFD Interior Optimization.nb` for 2nd-deriv A6 (and ditto for boundary). I want to optimize against what is really there, not what I infer from Kim07.

> RESPONSE: Here are a few lines of code which define the integrand for the spectral error function for first derivatives. If we define $\bar{\kappa}(\kappa) = \frac{N}{D}$ (numerator over denominator), then the first derivative integrand is $(\kappa D - N)^2$.
>
> `Num[\[Kappa]_] = 2*(a1*Sin[\[Kappa]] + a2*Sin[2*\[Kappa]] + a3*Sin[3*\[Kappa]]);`
>
> `Denom[\[Kappa]_] = 1 + 2*\[Alpha]*Cos[\[Kappa]] + 2*\[Beta]*Cos[2*\[Kappa]];`
>
> `\[Kappa]Bar[\[Kappa]_] = Num[\[Kappa]]/Denom[\[Kappa]];`
>
> `\[CapitalDelta]DH2[\[Kappa]_] = (\[Kappa]*Denom[\[Kappa]] - Num[\[Kappa]])^2 // Expand;`
>
> Here is the same information for second derivatives. If we define $\tilde{\kappa}(\kappa) = \frac{N}{D}$, then the second derivative integrand is $(\kappa^2 D - N)^2$.
>
> `Num[\[Kappa]_] = 4*(a1*(Sin[\[Kappa]/2])^2 + a2*(Sin[(2*\[Kappa])/2])^2 + a3*(Sin[(3*\[Kappa])/2])^2);`
>
> `Denom[\[Kappa]_] = 1 + 2*\[Alpha]*Cos[\[Kappa]] + 2*\[Beta]*Cos[2*\[Kappa]];`
>
> `\[Kappa]Bar[\[Kappa]_] = Num[\[Kappa]]/Denom[\[Kappa]];`
>
> `\[CapitalDelta]DH2[\[Kappa]_] = (\[Kappa]^2*Denom[\[Kappa]] - Num[\[Kappa]])^2 // Expand;`
>
> Note that all of the above was for the interior. For the boundaries, it is similar, but we use a complex norm since the spectral function $\bar{\kappa}(\kappa)$ or $\tilde{\kappa}(\kappa)$ is complex (includes dispersion and dissipation). For example, for second derivatives we have the integrand, which I denote by $\Delta^2$, is $\Delta^2 := (\kappa^2 - |\Re(\tilde{\kappa}(\kappa)) + i \Im(\tilde{\kappa}(\kappa))|)^2$. I don't want to copy and paste all the relevant code since it is a little bit messy, but I will point you to `CFD Boundary Optimization.nb`, then go to any operator Section (for example, under `Second Derivative`, go to the Section `P6 (A6) Scheme`), then go to the Subsection `Optimize $\hat{y}$`, then the Subsubsection `Node 2` (or any of the nodes), and the relevant code should be there. I can clear up remaining confusion if necessary.

---

## C. Topics I owe you a clearer explanation of in v1.2

I will fold these into the redraft. Listing here so nothing is forgotten.

### C1. The non-Dirichlet weighted boundary probe ($r=0.7$, Kim24 Eq. 17)

What it actually is, why Kim24 introduced it, whether it is genuinely problem-specific or whether it is a generic stability diagnostic. My preliminary view (see D1 below) is that it is the latter and may be exactly what catches your outlier eigenvalues — but I want to read Kim24 §3 carefully before defending that.

> RESPONSE: This would be helpful. I am still skeptical on your view in D1, but you are welcome to provide more justification for that claim or revise your claim. If possible, it may also be helpful to include ways in which we can test this claim. I also want to note that the non-Dirichlet parameter is introduced in Section 2.3, Eq. 17, so I believe Sections 2.3-3 provide good context into this parameter. After more carefully analyzing the Kim24 paper myself, I believe that implementing this $r$ parameter into our analysis of the eigenvalue spectra would not be too difficult; Kim24 mentions below Eq. 18 using the RHS matrix $\mathbf{Q}_r$ where the first column is replaced by $(1 - r) \mathbf{q}_0$, and we have essentially been using this approach with $r = 1$ already (setting the column to be all zeros), so again it doesn't seem difficult to implement. In Section 3, Kim24 mentions comparing $r = 1$ and $r = 0.7$ in Figure 2a. I see that using $r = 0.7$ essentially "flattens" the eigenvalue spectra more towards the imaginary axis which I could imagine coming closer to highlighting instabilities, but I don't think it will address the outlier eigenvalues. Still, I think let's address this topic in more depth in v1.2 and as mentioned before, let's think about ways to test this claim appropriately.

### C2. The role of $n$ and $\delta$ in Kim07's spectral functional

Why Kim included them, what mathematical role they play, and whether they could matter differently for $f''$ vs $f'$. I will derive this from Kim07 §2 rather than assert.

> RESPONSE: Yes, please address this. We had tested *just* the interior A4 scheme with/without the parameters $n$ and $\delta$ (the A4 scheme has identical structure to the Kim07 paper), and we saw minimal difference in the CFD coefficients when comparing with Table 2 in Kim07, which is why we discarded those parameters.

### C3. Whether the Zhou-Zhang weight shifts minima locations

I claimed in v1.1 A2 that the rational integrand is "cleaner." Your pushback is fair: a different weight will move the optimum unless the weight is a constant multiple of the original on the support. I will derive the stationary-point shift symbolically in v1.2 before recommending.

> RESPONSE: Yes, do this. I think this will help us to understand different options for the spectral error function, and that might give us a broader range of CFD operators that we can test as well.

### C4. Systematic vs arbitrary choice of integration limit $r$

Your $r \approx 0.85\pi$ was empirical. v1.2 will explain what $r$ controls (where on the modified-wavenumber curve we accept overshoot vs precision) and what determines a principled choice for $f''$ specifically.

> RESPONSE: Yes, do this. I should mention a recent test that we did after the v1.1 draft. For the A6 second derivative operators, I ran through the spectral error minimization procedure for the interior and boundary for $r \in \{0.60 \pi, 0.65 \pi, 0.70 \pi, 0.75 \pi, 0.80 \pi, 0.85 \pi\}$. Each of these gave a different set of about 40 operators (for example, 4 $\hat{y}$ minima in boundary node 2, 3 $\hat{y}$ minima in boundary node 1, and 4 $\hat{y}$ minima in boundary node 0, giving $4 \times 3 \times 4 = 48$ total operators to test). I arbitrarily chose 3 operators from each choice of the cutoff wavenumber $r$ in somewhat similar areas of the parameter space, and we tested each of these for accuracy and stability in an evolution. Out of the 18 total operators tested, 5 of them appeared to be stable (2 with $r = 0.60 \pi$, 1 with $r = 0.65 \pi$, and 2 with $r = 0.70 \pi$). Interestingly, they were stable in evolutions, but according to our eigenvalue stability criterion, these operators should not have been stable. Thus it seems that lowering the cutoff wavenumber as you mentioned was suggested in the Chen 2023 paper does help with stability. However, in our convergence tests, these A6 operators which were stable were converging at roughly 2nd order, which is not good. I would expect the accuracy to drop a little bit since we are optimizing the spectral error up to a lower wavenumber, but I would not expect the difference to be that dramatic. We will still be doing some more testing to understand this, but I think one thing I am interested in is how we effectively "translate" between a cutoff wavenumber $r \in (0, \pi)$ and something physical (for example, if we set $r = 0.65 \pi$ and looked at an electromagnetic multipole expansion, what is the highest order of multipole that we should expect to be well-resolved?).

### C5. CCD cost/benefit walkthrough

You want a clearer summary before deciding Track C. v1.2 will provide: precise per-step cost (twin-tridiagonal block elimination flop count), memory cost (2× per field × ~25 BSSN fields), and what specifically buys the ~8× $f''$ accuracy improvement (joint $f', f''$ unknowns vs decoupled).

> RESPONSE: Yes. We have started investigating the literature on the CCD method and it looks like it provides some genuine improvements worth pursuing, but we are still just in the beginning stages of understanding and testing this method. If you could explain it more carefully, synthesizing information from the Hixon 2017 paper, the Chu and Fan 1998 paper, and any other sources you find, that would be much appreciated.

### C6. Justifications (or removals) for each "skip" in v1.1 §5

The brief verdicts in v1.1 §5 were not enough. v1.2 will either give a real argument or remove the row.

> RESPONSE: Yes, provide more information as to why the ideas in Section 5 of v1.1 should not be pursued. Do *not* remove a row without asking.

---

## D. Things I think you may have missed (or where I want to push back)

### D1. The non-Dirichlet probe may not be problem-specific

Your v1.1 §6 Q1 response said you do not want to optimize to a specific test function, and that is why you reject the non-Dirichlet probe. My read is that the $r=0.7$ probe is a **stability analysis weighting**, not a test-function $L_\infty$ optimization — it is a way to expose eigenvalue behavior that pure-Dirichlet analysis hides. This may be exactly the diagnostic that catches your outlier eigenvalues. I want to defend or retract this in v1.2 after reading Kim24 §3 carefully.

> RESPONSE: I think you have misunderstood. Because we do not want to optimize to a specific test function, we are rejecting the $L_\infty$ optimization, not the non-Dirichlet probe (which is just a boundary condition). As for whether we should implement the non-Dirichlet probe, refer to Section C1 above for a more detailed description of how I think we ought to handle this.

### D2. Stability-during-optimization may be more tractable than you think

Brady-Livescu is brute-force ML, but **Deshpande's KKT framework** embeds stability as inequality constraints in a closed-form optimization. Mathematica's `NMinimize` with constraints does this natively. I would like to prototype it on the simplest case (e.g. one boundary row, fixed interior) before we conclude joint accuracy/stability optimization is implausible.

> RESPONSE: I think this is a reasonable pushback and is worth pursuing. Go ahead with prototyping a simple test case of this method. One thing I wanted to provide clarification for that I think fits best here; I mentioned in v1.1 that we have considered values $\hat{y} \in [0,5]$ for our boundary spectral error minimization. I would like to mention that in Eq. 24 of Kim24, Kim is considering values $\hat{y} \in (-6,6)$. Again, I have not seen many minima for negative values of $\hat{y}$, but I think this was worth mentioning.

### D3. Hybrid weighting may already give a clean answer for $f''$

If hybrids work for $f'$, the natural question is: what is the right hybrid weighting for $f''$ given that the modified-wavenumber target is $\kappa^2$, not $\kappa$? This is a symbolic derivation, not a numerical experiment, and I think it is worth doing **before** any wider-stencil or CCD work.

> RESPONSE: This point confuses me. We have verified that the hybrid method produces the anticipated qualitative changes in the spectral function for first derivatives and second derivatives. I don't understand what you mean by the "right hybrid weighting." Are you referring to the $\kappa$-dependence in the integrand of the spectral error function as mentioned in Section C3 above, or something else?

### D4. Verify Kim24's polynomial framework transfers to $f''$ unchanged

What differs between $f'$ and $f''$ in Kim24's boundary closure is the matched-derivative count and the row-constraint structure on the $f''$ side, not the polynomial framework. I want to confirm in v1.2 that your existing 2nd-deriv boundary notebook is using the framework this way (degree-13 polynomial, $n=8$ zero-derivative constraint, joint $\hat{y}$ optimization), vs. some other parameterization that drifted from Kim24.

> RESPONSE: In Kim24 Eq. 10, the constraints on the extrapolation polynomial for first derivatives are $P(j) = f_j$, $\frac{dP}{dx}(j) = f'_j$, and $\frac{d^n P}{dx^n}(\hat{y}) = 0$. In our approach for second derivatives, our conditions are $P(j) = f_j$, $\frac{d^2 P}{dx^2}(j) = f'_j$, and $\frac{d^n P}{dx^n}(\hat{y}) = 0$. We are still using a degree-13 polynomial. We have taken care to follow the advice in Kim24 below Eq. 10 that says $n$ should be larger than the desired order of accuracy of the system. We have varied the value of $n$, typically adopting a value of $n = 9$ for the A-schemes and a lower value such as $n = 6$ for the C-schemes. This is because the C-schemes only generate a degree-10 extrapolation polynomial due to having less parameters in the stencil. When we solve for the coefficients of the extrapolation polynomial, they are symbolic in terms of $\hat{y}$. What I have noticed is that the value of $n$ determines the dependence of the polynomial coefficients on $\hat{y}$. For example, if we take the A4-scheme as Kim24 does with $n = 8$, then the polynomial coefficients are quintic in $\hat{y}$ (contain $\hat{y}^5$ terms). What I have seen as a result is that for a lower value of $n$ (higher-order dependence on $\hat{y}$), the code takes a bit longer to run (because we are iteratively substituting in long expressions involving the polynomial coefficients $p_j$, so if the expressions for those coefficients are longer, the code takes a bit longer), but also a higher-order dependence on $\hat{y}$ gives more minima in the spectral error function for the boundary nodes. I have found it best to balance this by choosing a value of $n$ such that the polynomial coefficients are quartic in $\hat{y}$. When varying this, I have not seen qualitative differences (accuracy and stability) in the CFD operators which result, I have only seen a difference in the *number* of resulting operators we have to test. For example, I believe I did a test with first derivative A6 operators with $n \in \{10, 9, 8\}$ (cubic, quartic, and quintic in $\hat{y}$), found 144, 256, and 360 operators respectively, and *none* of them were stable according to the first derivative eigenvalue stability criterion. To summarize, the polynomial approach for second derivatives is the same as Kim24 but with the condition $\frac{d^2 P}{dx^2}(j) = f'_j$ replacing the condition $\frac{dP}{dx}(j) = f'_j$, we have tested varying values of $n$ in the $\hat{y}$ condition $\frac{d^n P}{dx^n}(\hat{y}) = 0$, and we have found that it does not change the accuracy or stability of our resulting operators.

### D5. The S6 (heptadiagonal LHS) attempt

You paused S6 because the eigenvalue spectra were unexpected. If 8th formal order is off the table for now (per your v1.1 §4 B2 response), S6 is also off the critical path. But: if the *eigenvalue surprises* on S6 share features with the *outlier eigenvalues* you have seen elsewhere, the diagnosis may be the same problem. Worth a brief look in v1.2 once B-scheme failure mode (B6) is characterized.

> RESPONSE: The S6 operators do not have outlier eigenvalues that were unexpected, the spectra just look similar to when I was working on the C6 operators, but I found a typo in that derivation that I had since corrected and which changed how the C6 operator eigenvalue spectra appeared visually. This caused me to believe that I had a typo in my derivation, which I will go through and check at some point. It is not related to Section B6 above. We can address this at a later point in the project, but for now let's not worry about the S6 operators.

### D6. Smoking-gun comparison still meaningful, just reframed

Even though my parameter-count claim in §2.1 was wrong, the Nate-vs-Tyler error gap **is** real and **is** evidence of the low-$\kappa$ vs total-error trade-off. v1.2 will keep the comparison but reframe it as "the trade-off the hybrids exist to address" rather than "lost free parameter."

> RESPONSE: This is correct. Refer to Section A6 above for more details on the hybrid approach.

---

## E. End-of-document open questions for you

### E1. Is there anything in v1.1 we have *not* discussed that you also want me to revisit?

(Sections, claims, references — anything you read and were unsure about but did not annotate.)

> RESPONSE: I believe it was partially addressed, but in v1.1 Section 2.4, I had asked the question, "What would cause the spectral error function for second derivatives to be more fragile than first derivatives?" Another thing we are interested at looking into is analyzing the spectral performance of these schemes in higher dimensions. For example, Figure 4 of Lele's paper shows polar plots of phase speed anisotropy, and so does Figure 1 of the Chu and Fan 1998 paper. I personally don't know a ton about these sorts of tests, but this is something we're interested in incorporating for our analysis. Also, due to switching which machine/account this project is running on, the directories to references have changed. Here is a summary of where to find the relevant literature we have discussed. In the directory `C:\Users\natha\Documents\School\Research\Hirschmann\ClaudeCFD\Papers`, I have the following papers, and I will try to refer to them by their names in this directory:
>
> - Chen 2023
> - Chu Fan CCD 1998
> - Deshpande 2019
> - Hixon CCD 2017
> - Jonathan Tyler Thesis (which we will also call the Tyler paper)
> - Kim07
> - Kim24
> - Lele 1992
> - Zhou Zhang 2011
>
> There are also a few more papers which we have not discussed but which I've included, mostly for my own reference:
>
> - Chu Fan 2000
> - Gurarslan CCD 2021
> - Zhang 2016
>
> Finally, I would also like you to search for additional relevant literature on the internet and provide it for our reference (but you do not have to read or analyze these additional references in-depth).

### E2. What is the success criterion for v1.2?

i.e. when v1.2 is ready, what should it let you decide? (Examples: "approve Track A and start coding", "decide whether to read the CCD papers in detail", "draft a paper outline".) Knowing the destination helps me write to it.

> RESPONSE: I would like v1.2 to address all items which were promised to be addressed in this document. It should provide a clear outline of which avenues/questions are left to address and pursue. For each of these, I would like a step-by-step of what to do next, including literature to read, derivations to follow/complete, test cases to construct and run, and results to analyze and interpret.

### E3. Timeframe / cadence

Are we drafting v1.2 in one pass after you finish reviewing this file, or iteratively (you respond to a few items, I draft partial sections, repeat)?

> RESPONSE: I have responded to every item in this follow-up, so please draft v1.2 all at once. If there are additional things that need review before doing so, let me know. Please ask if anything remains unclear or if more information is needed.
