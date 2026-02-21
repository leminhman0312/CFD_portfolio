import numpy as np


# ============================================================================
# Thomas tridiagonal solver
# This matches the C++ function signature and indexing style
# Solves for u[1..N] given a,b,c,d arrays indexed 1..N
# Here:
#   dprime[i] = d[i] - (b[i]*a[i-1])/dprime[i-1]
#   cprime[i] = c[i] - (cprime[i-1]*b[i])/dprime[i-1]
# ============================================================================
def thomasTriDiagonal(imax: int, a: np.ndarray, b: np.ndarray, c: np.ndarray,
                      d: np.ndarray, u: np.ndarray) -> None:
    dprime = np.zeros(imax + 1, dtype=float)
    cprime = np.zeros(imax + 1, dtype=float)

    dprime[1] = d[1]
    cprime[1] = c[1]

    for i in range(2, imax + 1):
        dprime[i] = d[i] - (b[i] * a[i - 1]) / dprime[i - 1]
        cprime[i] = c[i] - (cprime[i - 1] * b[i]) / dprime[i - 1]

    u[imax] = cprime[imax] / dprime[imax]

    for i in range(imax - 1, 0, -1):
        u[i] = (cprime[i] - a[i] * u[i + 1]) / dprime[i]
