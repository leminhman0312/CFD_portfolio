import numpy as np
from linalg import thomasTriDiagonal


# ============================================================================
# Numerical schemes
# All arrays use length imax+1 and indices 1..imax are valid
# ============================================================================
def FTCS_explicit(nmax: int, F_num: float, tboundary: float, t0: float,
                  imax: int, out_un: np.ndarray) -> None:
    u = np.zeros(imax + 1, dtype=float)
    un = np.zeros(imax + 1, dtype=float)

    for i in range(2, imax):
        u[i] = t0
        un[i] = t0

    u[1] = tboundary
    u[imax] = tboundary
    un[1] = tboundary
    un[imax] = tboundary

    for _n in range(1, nmax + 1):
        for i in range(1, imax + 1):
            u[i] = un[i]

        for i in range(2, imax):
            un[i] = u[i] + F_num * (u[i + 1] - 2.0 * u[i] + u[i - 1])

        un[1] = tboundary
        un[imax] = tboundary

    for i in range(1, imax + 1):
        out_un[i] = un[i]


def FTCS_return_implicit_interior(nmax: int, F_num: float, tboundary: float,
                                  t0: float, imax: int, u_full: np.ndarray) -> None:
    N = imax - 2  # interior unknown count

    a = np.zeros(N + 1, dtype=float)
    b = np.zeros(N + 1, dtype=float)
    c = np.zeros(N + 1, dtype=float)
    d = np.zeros(N + 1, dtype=float)
    u_int = np.zeros(N + 1, dtype=float)

    u_full[1] = tboundary
    u_full[imax] = tboundary
    for i in range(2, imax):
        u_full[i] = t0

    for j in range(1, N + 1):
        d[j] = 1.0 + 2.0 * F_num
        a[j] = -F_num
        b[j] = -F_num

    b[1] = 0.0
    a[N] = 0.0

    for _t in range(1, nmax + 1):
        for j in range(1, N + 1):
            i = j + 1
            c[j] = u_full[i]

        c[1] += F_num * tboundary
        c[N] += F_num * tboundary

        thomasTriDiagonal(N, a, b, c, d, u_int)

        u_full[1] = tboundary
        u_full[imax] = tboundary
        for j in range(1, N + 1):
            i = j + 1
            u_full[i] = u_int[j]


def DufortFrankel(nmax: int, F_num: float, tboundary: float, t0: float,
                  imax: int, out_un: np.ndarray) -> None:
    dval = 2.0 * F_num

    un_m1 = np.zeros(imax + 1, dtype=float)
    un_p1 = np.zeros(imax + 1, dtype=float)
    un = np.zeros(imax + 1, dtype=float)

    un_m1[1] = tboundary
    un_m1[imax] = tboundary
    for i in range(2, imax):
        un_m1[i] = t0

    u1 = np.zeros(imax + 1, dtype=float)
    FTCS_return_implicit_interior(1, F_num, tboundary, t0, imax, u1)
    for i in range(1, imax + 1):
        un[i] = u1[i]

    for _n in range(1, nmax):
        for i in range(2, imax):
            un_p1[i] = ((1.0 - dval) * un_m1[i] + dval * (un[i + 1] + un[i - 1])) / (1.0 + dval)

        un_p1[1] = tboundary
        un_p1[imax] = tboundary

        for i in range(1, imax + 1):
            un_m1[i] = un[i]
            un[i] = un_p1[i]

    for i in range(1, imax + 1):
        out_un[i] = un[i]


def FTCS_implicit(nmax: int, F_num: float, tboundary: float, t0: float,
                  imax: int, out_un: np.ndarray) -> None:
    a = np.zeros(imax + 1, dtype=float)
    b = np.zeros(imax + 1, dtype=float)
    c = np.zeros(imax + 1, dtype=float)
    d = np.zeros(imax + 1, dtype=float)
    un = np.zeros(imax + 1, dtype=float)

    un[1] = tboundary
    un[imax] = tboundary
    for i in range(2, imax):
        un[i] = t0

    d[1] = 1.0
    a[1] = 0.0
    b[1] = 0.0

    d[imax] = 1.0
    a[imax] = 0.0
    b[imax] = 0.0

    for i in range(2, imax):
        d[i] = 1.0 + 2.0 * F_num
        a[i] = -F_num
        b[i] = -F_num

    for _t in range(1, nmax + 1):
        c[1] = tboundary
        c[imax] = tboundary
        for i in range(2, imax):
            c[i] = un[i]

        thomasTriDiagonal(imax, a, b, c, d, un)

    for i in range(1, imax + 1):
        out_un[i] = un[i]


def CrankNicolson(nmax: int, F_num: float, tboundary: float, t0: float,
                  imax: int, out_un: np.ndarray) -> None:
    d_NC = F_num / 2.0

    a = np.zeros(imax + 1, dtype=float)
    b = np.zeros(imax + 1, dtype=float)
    c = np.zeros(imax + 1, dtype=float)
    d = np.zeros(imax + 1, dtype=float)

    u0 = np.zeros(imax + 1, dtype=float)
    u_half = np.zeros(imax + 1, dtype=float)

    u0[1] = tboundary
    u0[imax] = tboundary
    for i in range(2, imax):
        u0[i] = t0

    u_half[1] = tboundary
    u_half[imax] = tboundary

    d[1] = 1.0
    a[1] = 0.0
    b[1] = 0.0

    d[imax] = 1.0
    a[imax] = 0.0
    b[imax] = 0.0

    for i in range(2, imax):
        d[i] = 1.0 + 2.0 * d_NC
        a[i] = -d_NC
        b[i] = -d_NC

    for _n in range(1, nmax + 1):
        u_half[1] = tboundary
        u_half[imax] = tboundary

        for i in range(2, imax):
            u_half[i] = u0[i] + d_NC * (u0[i + 1] - 2.0 * u0[i] + u0[i - 1])

        c[1] = tboundary
        c[imax] = tboundary
        for i in range(2, imax):
            c[i] = u_half[i]

        thomasTriDiagonal(imax, a, b, c, d, u0)

        u0[1] = tboundary
        u0[imax] = tboundary

    for i in range(1, imax + 1):
        out_un[i] = u0[i]
