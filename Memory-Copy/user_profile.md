---
name: User profile — Nate Garey
description: Nate's role, background, and working style for the ClaudeCFD project
type: user
---

Nate Garey (ngarey@byu.edu) is a graduate student at Brigham Young University (BYU) working in the numerical relativity group under Prof. Hirschmann. His primary research project is developing optimized compact finite difference (CFD) operators for use in the Dendro-GR BSSN code.

**Technical background:** Strong in Mathematica and numerical methods; familiar with compact FD literature (Kim07, Kim24, Lele 1992, Tyler thesis). Working knowledge of numerical relativity (BSSN formalism, gravitational waves). Less familiar with newer CCD literature (Chu-Fan, Hixon) and some stability analysis details (non-Dirichlet probe).

**Role in the project:** Nate implements all CFD derivations and Mathematica notebooks himself. His collaborator Jonathan Tyler (BYU Master's thesis) developed an alternative spectral optimization approach ("Tyler method"). A previous Claude Code session on a different machine/account produced v1.1 of the analysis; Nate is continuing it on this machine.

**Preferences:** Prefers spectral error optimization over test-function-dependent approaches (e.g., L∞ on exp(sin x)). Wants operators that are general-purpose across Maxwell, NSM, Teukolsky, and eventually BSSN — not tuned to a specific test function. Skeptical of problem-specific stability diagnostics until clearly justified.
