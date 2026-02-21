import math
import numpy as np


# ============================================================================
# Exact solution (series), matches the C++ exact_solution_profile
# x[1] = 0, x[i] = (i-1)*dx, L = (imax-1)*dx
# ============================================================================
def exact_solution_profile(imax: int, x: np.ndarray, t: float,
                           alpha: float, L: float,
                           tboundary: float, t0: float,
                           nterms: int, out_T: np.ndarray) -> None:
    pi = 4.0 * math.atan(1.0)
    A = (t0 - tboundary)

    for i in range(1, imax + 1):
        s = 0.0
        xi = x[i]
        for n in range(1, nterms + 1):
            coeff = (2.0 * A / (n * pi)) * (1.0 - ((-1.0) ** n))
            if coeff == 0.0:
                continue
            k = n * pi / L
            s += coeff * math.sin(k * xi) * math.exp(-alpha * k * k * t)
        out_T[i] = tboundary + s

    out_T[1] = tboundary
    out_T[imax] = tboundary
