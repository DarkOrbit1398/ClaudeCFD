---
name: ClaudeCFD project state
description: Current state of the compact finite difference operator development project for Dendro-GR
type: project
---

**Goal:** Develop accurate and stable 6th-order compact finite difference (CFD) operators for second derivatives, for use in Dendro-GR (octree-based AMR NR code evolving BSSN with 24 fields, 13×13×13 blocks, k=3 ghost zones).

**Key context:**
- Interior optimization follows Kim07 (spectral error minimization via ∂E/∂α=0); we discard Kim07's n and δ parameters (tested: negligible effect).
- Boundary closure follows Kim24 (degree-13 extrapolating polynomial, per-row ŷ values, n=9 for 6th-order A-schemes).
- We use spectral error optimization for both interior AND boundary (not Kim24's L∞ for boundary) — intentional, for test-function independence.
- Current best: A6 2nd-derivative operators exist but converge at explicit-FD-equivalent accuracy. Root cause is the low-κ vs total-κ trade-off (Kim-like vs Tyler-like optimization), not a free-parameter deficit.
- Hybrid schemes (mixing Kim-like ∂E/∂α=0 and Tyler-like ω̃(ωᵢ)=ωᵢ constraints) are implemented for 1st derivatives and are a first-class option for 2nd derivatives.

**Notebook paths:** `C:\Users\natha\Documents\School\Research\Hirschmann\HirschmannGroupBYU\Notes\MathematicaNotebooks\`
- `CFD Interior Optimization.nb`, `CFD Boundary Optimization.nb`, `CFD Eigenvalue Stability.nb`

**Papers path:** `C:\Users\natha\Documents\School\Research\Hirschmann\ClaudeCFD\Papers\`

**Planning documents:** `C:\Users\natha\Documents\School\Research\Hirschmann\ClaudeCFD\Plans\`
- `project.md` (original brief), `analysis_v1.1.md` (first analysis, annotated by Nate), `review_followup.md` (Nate's responses to v1.1), `analysis_v1.2.md` (current analysis, 2026-05-07)

**Current notebook issues (as of 2026-05-07):**
- B4/C4/B6/C6 2nd-deriv work and hybrid interior schemes may be missing due to Git overwrite — to be restored.
- `yHat1List[[5]]` index bug in 2nd-deriv P6(A6) — placeholder fix needed.
- S6 (heptadiagonal LHS) paused, suspected derivation typo.

**Active tracks (v1.2):**
- Track A: Tighten existing A6 2nd-deriv (r scan, non-Dirichlet probe, NMinimize for ŷ, hybrid interior, Zhou-Zhang weight). **Priority diagnosis: why do stable operators at r≤0.70π converge at ~2nd order?**
- Track B: Diagnose B-scheme failure (low convergence order despite good spectral resolution).
- Track C: CCD approach (Chu-Fan 1998 + Hixon 2017). **Must read Daszuta et al. 2024 before starting.**

**Key facts:**
- Stability criterion for 2nd deriv: Re(μ) ≤ 0, Im(μ) = 0 for eigenvalues μ of P⁻¹Q; equivalently Re(λ)=0, Im(λ)≥0 for λ=√μ.
- Currently NOT using Kim24 non-Dirichlet probe (r=0.7); implementation is easy and is planned in Track A.
- ŷ scan range: currently [0,5] in steps 0.03; should extend to (-6,6) per Kim24.
- Dendro-GR has 24 BSSN fields (not 25), 13×13×13 blocks confirmed.
- Brady & Livescu 2019 is NOT an ML paper — uses direct nonlinear testing on Euler equations.
- Deshpande 2021 (preprint 2019): KKT closed-form optimization for periodic interior schemes.
- Daszuta et al. 2024 (*J. Comput. Phys.*): first NR compact FD paper — must obtain and read.
