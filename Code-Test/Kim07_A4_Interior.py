"""
Kim07 interior A4 optimization via KKT linear system.
Reproduces Tables 1 and 2 of Kim (2007), AIAA J. 45(1).

Scheme (1st derivative, interior):
  β f'_{i-2} + α f'_{i-1} + f'_i + α f'_{i+1} + β f'_{i+2}
      = (1/h) [a1(f_{i+1}-f_{i-1}) + a2(f_{i+2}-f_{i-2}) + a3(f_{i+3}-f_{i-3})]

4th-order Taylor constraints (matching κ¹ and κ³ in N = κD near κ=0):
  -2α - 2β + 2a1 + 4a2 + 6a3 = 1   (Kim07 Eq 2)
  -3α - 12β + a1 + 8a2 + 27a3 = 0  (Kim07 Eq 3)

Spectral error functional (Kim07 Eq 7):
  E = ∫₀ʳ [(1+δ)κ D(κ) - N(κ)]² (κ/r)ⁿ dκ
  D(κ) = 1 + 2α cos κ + 2β cos 2κ
  N(κ) = 2[a1 sin κ + a2 sin 2κ + a3 sin 3κ]

E is quadratic in x = [α, β, a1, a2, a3] and constraints are linear, so the
constrained minimum solves the 7×7 KKT system [[M, Cᵀ],[C, 0]][x;λ] = [-f; q]
where M[i,j] = ∫φᵢφⱼ w dκ, f[i] = ∫φ₀φᵢ w dκ, w = (κ/r)ⁿ.

Basis functions (coefficients in the residual (1+δ)κD - N):
  φ₀    = (1+δ)κ
  φ_α   = 2(1+δ)κ cos κ
  φ_β   = 2(1+δ)κ cos 2κ
  φ_a1  = -2 sin κ
  φ_a2  = -2 sin 2κ
  φ_a3  = -2 sin 3κ
"""

import numpy as np
from scipy.integrate import quad


def _phi(k, d):
    """Basis vector [phi_alpha, phi_beta, phi_a1, phi_a2, phi_a3] at wavenumber k."""
    c = 1.0 + d
    return np.array([
         2.0 * c * k * np.cos(k),
         2.0 * c * k * np.cos(2.0 * k),
        -2.0 * np.sin(k),
        -2.0 * np.sin(2.0 * k),
        -2.0 * np.sin(3.0 * k),
    ])


def solve_kim07_a4(r, delta, n):
    """
    Solve Kim07 A4 interior optimization for given (r, delta, n).

    Parameters
    ----------
    r     : upper wavenumber limit for integration (spectral resolution parameter)
    delta : small spectral offset in the functional target (Kim07 Eq 7)
    n     : weight exponent (Kim07 Eq 7)

    Returns
    -------
    dict with keys: r, delta, n, alpha, beta, a1, a2, a3
    """
    d = delta
    phi0_val = lambda k: (1.0 + d) * k

    def w(k):
        return (k / r) ** n

    # Gram matrix  M[i,j] = ∫₀ʳ phi_i(k) phi_j(k) w(k) dk
    M = np.zeros((5, 5))
    for i in range(5):
        for j in range(i, 5):
            val, _ = quad(
                lambda k, i=i, j=j: _phi(k, d)[i] * _phi(k, d)[j] * w(k),
                0.0, r, limit=200,
            )
            M[i, j] = M[j, i] = val

    # RHS vector  f[i] = ∫₀ʳ phi0(k) phi_i(k) w(k) dk
    f = np.array([
        quad(lambda k, i=i: phi0_val(k) * _phi(k, d)[i] * w(k),
             0.0, r, limit=200)[0]
        for i in range(5)
    ])

    # Taylor constraint matrix  C x = q,  x = [α, β, a1, a2, a3]
    C = np.array([
        [-2.0, -2.0,  2.0,  4.0,  6.0],
        [-3.0, -12.0, 1.0,  8.0, 27.0],
    ])
    q = np.array([1.0, 0.0])

    # KKT system  [[M, Cᵀ], [C, 0]] · [x; λ] = [-f; q]
    K = np.block([[M, C.T], [C, np.zeros((2, 2))]])
    rhs = np.concatenate([-f, q])
    sol = np.linalg.solve(K, rhs)

    alpha, beta, a1, a2, a3 = sol[:5]
    return dict(r=r, delta=delta, n=n,
                alpha=alpha, beta=beta, a1=a1, a2=a2, a3=a3)


def resolution_limit(coeff, tol=1e-3, n_pts=10000):
    """Return κ_c / π where |κ̄/κ - 1| first exceeds tol (default 0.1%)."""
    α, β = coeff['alpha'], coeff['beta']
    a1, a2, a3 = coeff['a1'], coeff['a2'], coeff['a3']
    k = np.linspace(1e-8, np.pi, n_pts)
    D = 1.0 + 2.0*α*np.cos(k) + 2.0*β*np.cos(2.0*k)
    N = 2.0*a1*np.sin(k) + 2.0*a2*np.sin(2.0*k) + 2.0*a3*np.sin(3.0*k)
    exceed = np.where(np.abs(N / D / k - 1.0) > tol)[0]
    return k[exceed[0]] / np.pi if len(exceed) else 1.0


def print_result(coeff, label=''):
    if label:
        print(f'\n--- {label} ---')
    print(f"  alpha = {coeff['alpha']:.16f}")
    print(f"  beta  = {coeff['beta']:.16f}")
    print(f"  a1    = {coeff['a1']:.16f}")
    print(f"  a2    = {coeff['a2']:.16f}")
    print(f"  a3    = {coeff['a3']:.16f}")
    print(f"  kappa_c^0.1% = {resolution_limit(coeff):.4f} pi")


# ---------------------------------------------------------------------------
# Reference values from Kim07 Table 1
# ---------------------------------------------------------------------------
TABLE1_REF = dict(
    alpha = 0.5862704032801503,
    beta  = 0.0954953355017055,
    a1    = 0.6431406736919156,
    a2    = 0.2586011023495066,
    a3    = 0.007140953479797375,
)

# ---------------------------------------------------------------------------
# Cases from Kim07 Table 2  (r, delta, n, label)
# ---------------------------------------------------------------------------
TABLE2_CASES = [
    dict(r=2.661, delta=-0.000210, n=15, label='r=2.661  delta=-0.000210  n=15'),
    dict(r=2.671, delta=-0.000246, n=11, label='r=2.671  delta=-0.000246  n=11'),
    dict(r=2.669, delta=-0.000260, n=10, label='r=2.669  delta=-0.000260  n=10'),
    dict(r=2.672, delta=-0.000233, n=10, label='r=2.672  delta=-0.000233  n=10'),
]


if __name__ == '__main__':
    print('=' * 60)
    print('Kim07 Interior A4 Optimization  (KKT linear system)')
    print('=' * 60)

    # ------------------------------------------------------------------
    print('\n=== TABLE 1: r=2.672, delta=-0.000233, n=10 ===')
    res1 = solve_kim07_a4(r=2.672, delta=-0.000233, n=10)
    print_result(res1, 'computed')

    print('\nKim07 Table 1 reference:')
    for k, v in TABLE1_REF.items():
        print(f'  {k} = {v:.16f}')

    print('\nErrors vs Kim07 Table 1:')
    for k in ('alpha', 'beta', 'a1', 'a2', 'a3'):
        print(f'  |{k} - ref| = {abs(res1[k] - TABLE1_REF[k]):.2e}')

    # ------------------------------------------------------------------
    print('\n=== TABLE 2 ===')
    for case in TABLE2_CASES:
        res = solve_kim07_a4(r=case['r'], delta=case['delta'], n=case['n'])
        print_result(res, case['label'])
