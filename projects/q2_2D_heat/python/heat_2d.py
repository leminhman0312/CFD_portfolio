#!/usr/bin/env python3
"""
2D Heat Conduction Solver in Python

Solves the transient two dimensional heat equation

    ∂T/∂t = α ( ∂²T/∂x² + ∂²T/∂y² )

Numerical methods
    FTCS explicit scheme
    Implicit ADI scheme using Thomas tridiagonal solvers

Run modes
    python heat2d.py sim
    python heat2d.py convergence
    python heat2d.py all

Directories created automatically
    data
    plot
    gnuplot_scripts
"""

import math
import os
import sys
import subprocess


# ----------------------------------------------------------------------------
# Directory helpers
# ----------------------------------------------------------------------------

def ensure_dir(name: str) -> None:
    os.makedirs(name, exist_ok=True)


def ensure_project_dirs() -> None:
    ensure_dir("data")
    ensure_dir("plot")
    ensure_dir("gnuplot_scripts")


# ----------------------------------------------------------------------------
# Error norms
# ----------------------------------------------------------------------------

def error_L2(u, uref) -> float:
    imax = len(u) - 1
    jmax = len(u[1]) - 1

    total = 0.0
    n = 0

    for i in range(2, imax):
        for j in range(2, jmax):
            e = u[i][j] - uref[i][j]
            total += e * e
            n += 1

    return math.sqrt(total / n) if n > 0 else 0.0


def error_Linf(u, uref) -> float:
    imax = len(u) - 1
    jmax = len(u[1]) - 1

    max_err = 0.0

    for i in range(2, imax):
        for j in range(2, jmax):
            e = abs(u[i][j] - uref[i][j])
            if e > max_err:
                max_err = e

    return max_err


# ----------------------------------------------------------------------------
# High level simulation driver
# ----------------------------------------------------------------------------

def sim(u0,
        t_end,
        deltax, deltay,
        alpha,
        dt_implicit,
        dt_explicit_given,
        t1, t2, t3, t4) -> None:

    print("\nSIMULATING")

    write_field_xyz("data/initial.dat", u0, deltax, deltay)
    plotContourMatlabLike(
        "data/initial.dat",
        "plot/contour_initial.png",
        0.0,
        "Initial conditions"
    )

    nmax_implicit = round(t_end / dt_implicit)

    u_implicit = FTCS_implicit_ADI(
        u0, nmax_implicit,
        deltax, deltay, dt_implicit,
        alpha, t1, t2, t3, t4
    )

    write_field_xyz("data/implicit.dat", u_implicit, deltax, deltay)
    plotContourMatlabLike(
        "data/implicit.dat",
        "plot/contour_implicit.png",
        t_end,
        "Implicit ADI"
    )

    fx = alpha * dt_explicit_given / (deltax * deltax)
    fy = alpha * dt_explicit_given / (deltay * deltay)
    s = fx + fy

    if s <= 0.5:
        u_explicit = FTCS_Explicit(
            u0, t_end, deltax, deltay,
            dt_explicit_given, alpha,
            t1, t2, t3, t4
        )

        write_field_xyz("data/explicit.dat", u_explicit, deltax, deltay)
        plotContourMatlabLike(
            "data/explicit.dat",
            "plot/contour_explicit.png",
            t_end,
            "Explicit FTCS"
        )
    else:
        dt_safe = 0.5 / (alpha * (1.0 / (deltax * deltax) + 1.0 / (deltay * deltay)))

        u_fail = FTCS_Explicit(
            u0, t_end, deltax, deltay,
            dt_explicit_given, alpha,
            t1, t2, t3, t4
        )

        u_safe = FTCS_Explicit(
            u0, t_end, deltax, deltay,
            dt_safe, alpha,
            t1, t2, t3, t4
        )

        write_field_xyz("data/explicit_failed.dat", u_fail, deltax, deltay)
        write_field_xyz("data/explicit_safe.dat", u_safe, deltax, deltay)

        plotContourMatlabLike(
            "data/explicit_failed.dat",
            "plot/contour_explicit_failed.png",
            t_end,
            "Explicit unstable dt"
        )

        plotContourMatlabLike(
            "data/explicit_safe.dat",
            "plot/contour_explicit_safe.png",
            t_end,
            "Explicit safe dt"
        )


# ----------------------------------------------------------------------------
# Time convergence study
# ----------------------------------------------------------------------------

def convergence(u0,
                t_end,
                deltax, deltay,
                alpha,
                t1, t2, t3, t4) -> None:

    print("\nTIME CONVERGENCE STUDY (implicit ADI)")
    print("dt        L2 error      Linf error")
    print("----------------------------------")

    dt_start = 0.04
    dt_ratio = 0.5
    dt_min = 1.0e-6

    dt_list = []
    dt = dt_start
    while dt >= dt_min:
        dt_list.append(dt)
        dt *= dt_ratio

    dt_ref = dt_list[-1]
    nmax_ref = round(t_end / dt_ref)

    u_ref = FTCS_implicit_ADI(
        u0, nmax_ref,
        deltax, deltay, dt_ref,
        alpha, t1, t2, t3, t4
    )

    with open("data/time_convergence.dat", "w") as f:
        for dt in dt_list:
            nmax = round(t_end / dt)

            u = FTCS_implicit_ADI(
                u0, nmax,
                deltax, deltay, dt,
                alpha, t1, t2, t3, t4
            )

            L2 = error_L2(u, u_ref)
            Linf = error_Linf(u, u_ref)

            print(f"{dt:<8.5f}  {L2:.6e}  {Linf:.6e}")
            f.write(f"{dt:.10f} {L2:.10e} {Linf:.10e}\n")

    print("\nWrote data/time_convergence.dat")
    plot_time_convergence("data/time_convergence.dat", "plot/time_convergence.png")


# ----------------------------------------------------------------------------
# Initialize field
# ----------------------------------------------------------------------------

def initializeField(imax, jmax, t0, t1, t2, t3, t4):
    """
    Uses 1 based indexing style to match the C++ code.
    Index 0 is unused.
    """
    u = [[t0 for _ in range(jmax + 1)] for _ in range(imax + 1)]

    for j in range(1, jmax + 1):
        u[1][j] = t2
        u[imax][j] = t4

    for i in range(1, imax + 1):
        u[i][1] = t1
        u[i][jmax] = t3

    return u


# ----------------------------------------------------------------------------
# Explicit FTCS solver
# ----------------------------------------------------------------------------

def FTCS_Explicit(u0,
                  t_end,
                  deltax, deltay,
                  dt,
                  alpha,
                  t1, t2, t3, t4):

    u = [row[:] for row in u0]

    imax = len(u) - 1
    jmax = len(u[1]) - 1

    u_new = [[0.0 for _ in range(jmax + 1)] for _ in range(imax + 1)]

    nfull = math.floor(t_end / dt)
    dt_last = t_end - nfull * dt
    nsteps = nfull + 1 if dt_last > 0.0 else nfull

    for n in range(1, nsteps + 1):
        dt_n = dt_last if (n == nsteps and dt_last > 0.0) else dt

        fx = alpha * dt_n / (deltax * deltax)
        fy = alpha * dt_n / (deltay * deltay)

        for i in range(1, imax + 1):
            for j in range(1, jmax + 1):
                u_new[i][j] = u[i][j]

        for i in range(2, imax):
            for j in range(2, jmax):
                u_new[i][j] = (
                    (1.0 - 2.0 * fx - 2.0 * fy) * u[i][j]
                    + fx * (u[i + 1][j] + u[i - 1][j])
                    + fy * (u[i][j + 1] + u[i][j - 1])
                )

        for j in range(1, jmax + 1):
            u_new[1][j] = t2
            u_new[imax][j] = t4

        for i in range(1, imax + 1):
            u_new[i][1] = t1
            u_new[i][jmax] = t3

        u, u_new = u_new, u

    return u


# ----------------------------------------------------------------------------
# Implicit ADI solver
# ----------------------------------------------------------------------------

def FTCS_implicit_ADI(u0,
                      nmax,
                      deltax, deltay,
                      dt,
                      alpha,
                      t1, t2, t3, t4):

    u = [row[:] for row in u0]

    imax = len(u) - 1
    jmax = len(u[1]) - 1

    u_dummy = [[0.0 for _ in range(jmax + 1)] for _ in range(imax + 1)]

    fx = alpha * dt / (deltax * deltax)
    fy = alpha * dt / (deltay * deltay)

    ax = [0.0] * (imax + 1)
    bx = [0.0] * (imax + 1)
    cx = [0.0] * (imax + 1)
    dx = [0.0] * (imax + 1)

    for i in range(2, imax):
        ax[i] = -fx / 2.0
        bx[i] = -fx / 2.0
        dx[i] = 1.0 + fx
    dx[1] = 1.0
    dx[imax] = 1.0

    ay = [0.0] * (jmax + 1)
    by = [0.0] * (jmax + 1)
    cy = [0.0] * (jmax + 1)
    dy = [0.0] * (jmax + 1)

    for j in range(2, jmax):
        ay[j] = -fy / 2.0
        by[j] = -fy / 2.0
        dy[j] = 1.0 + fy
    dy[1] = 1.0
    dy[jmax] = 1.0

    for _ in range(1, nmax + 1):
        for j in range(1, jmax + 1):
            u[1][j] = t2
            u[imax][j] = t4

        for i in range(1, imax + 1):
            u[i][1] = t1
            u[i][jmax] = t3

        # X sweep
        for j in range(2, jmax):
            for i in range(2, imax):
                cx[i] = (
                    (1.0 - fy) * u[i][j]
                    + (fy / 2.0) * (u[i][j + 1] + u[i][j - 1])
                )

            cx[2] += (fx / 2.0) * u[1][j]
            cx[imax - 1] += (fx / 2.0) * u[imax][j]

            solx = [0.0] * (imax + 1)
            solx[1] = u[1][j]
            solx[imax] = u[imax][j]

            cx[1] = u[1][j]
            cx[imax] = u[imax][j]

            thomasTriDiagonal(imax, ax, bx, cx, dx, solx)

            for i in range(1, imax + 1):
                u_dummy[i][j] = solx[i]

        for j in range(1, jmax + 1):
            u_dummy[1][j] = t2
            u_dummy[imax][j] = t4

        for i in range(1, imax + 1):
            u_dummy[i][1] = t1
            u_dummy[i][jmax] = t3

        # Y sweep
        for i in range(2, imax):
            for j in range(2, jmax):
                cy[j] = (
                    (1.0 - fx) * u_dummy[i][j]
                    + (fx / 2.0) * (u_dummy[i + 1][j] + u_dummy[i - 1][j])
                )

            cy[2] += (fy / 2.0) * u_dummy[i][1]
            cy[jmax - 1] += (fy / 2.0) * u_dummy[i][jmax]

            soly = [0.0] * (jmax + 1)
            soly[1] = u_dummy[i][1]
            soly[jmax] = u_dummy[i][jmax]

            cy[1] = u_dummy[i][1]
            cy[jmax] = u_dummy[i][jmax]

            thomasTriDiagonal(jmax, ay, by, cy, dy, soly)

            for j in range(1, jmax + 1):
                u[i][j] = soly[j]

        for j in range(1, jmax + 1):
            u[1][j] = t2
            u[imax][j] = t4

        for i in range(1, imax + 1):
            u[i][1] = t1
            u[i][jmax] = t3

    return u


# ----------------------------------------------------------------------------
# Thomas tridiagonal solver
# ----------------------------------------------------------------------------

def thomasTriDiagonal(n, a, b, c, d, u):
    dprime = [0.0] * (n + 1)
    cprime = [0.0] * (n + 1)

    dprime[1] = d[1]
    cprime[1] = c[1]

    for i in range(2, n + 1):
        dprime[i] = d[i] - (b[i] * a[i - 1]) / dprime[i - 1]
        cprime[i] = c[i] - (cprime[i - 1] * b[i]) / dprime[i - 1]

    u[n] = cprime[n] / dprime[n]

    for i in range(n - 1, 0, -1):
        u[i] = (cprime[i] - a[i] * u[i + 1]) / dprime[i]


# ----------------------------------------------------------------------------
# Write field to file
# ----------------------------------------------------------------------------

def write_field_xyz(filename, u, deltax, deltay) -> None:
    imax = len(u) - 1
    jmax = len(u[1]) - 1

    with open(filename, "w") as f:
        for j in range(1, jmax + 1):
            y = (j - 1) * deltay
            for i in range(1, imax + 1):
                x = (i - 1) * deltax
                f.write(f"{x:.8f} {y:.8f} {u[i][j]:.10f}\n")
            f.write("\n")


# ----------------------------------------------------------------------------
# Plot helpers
# ----------------------------------------------------------------------------

def plotContourMatlabLike(datafile, outpng, time_hr, scheme) -> None:
    cmd = [
        "gnuplot",
        "-e",
        f'datafile="{datafile}"; outpng="{outpng}"; tlabel={time_hr:.3f}; scheme="{scheme}"',
        "gnuplot_scripts/plot_contour_2d_matlab_like.gp"
    ]

    print(f"\n{scheme}")
    try:
        subprocess.run(cmd, check=True)
    except FileNotFoundError:
        print("gnuplot not found. Skipping contour plot.")
    except subprocess.CalledProcessError as e:
        print(f"gnuplot failed: {e}")


def plot_time_convergence(datafile, outpng) -> None:
    cmd = [
        "gnuplot",
        "-e",
        f'datafile="{datafile}"; outpng="{outpng}"',
        "gnuplot_scripts/plot_convergence.gp"
    ]

    print("\nTime convergence plot")
    try:
        subprocess.run(cmd, check=True)
    except FileNotFoundError:
        print("gnuplot not found. Skipping convergence plot.")
    except subprocess.CalledProcessError as e:
        print(f"gnuplot failed: {e}")


# ----------------------------------------------------------------------------
# Main
# ----------------------------------------------------------------------------

def main():
    ensure_project_dirs()

    arg1 = sys.argv[1] if len(sys.argv) > 1 else None

    do_sim = (len(sys.argv) == 1) or (arg1 == "sim")
    do_conv = (arg1 == "convergence")
    do_all = (arg1 == "all")

    t_end = 0.5

    deltax = 0.1
    deltay = 0.1
    alpha = 0.645

    dt_implicit = 0.01
    dt_explicit_given = 0.01

    xmin, xmax = 0.0, 3.5
    ymin, ymax = 0.0, 3.5

    imax = math.ceil((xmax - xmin) / deltax + 1.0)
    jmax = math.ceil((ymax - ymin) / deltay + 1.0)

    t0 = 0.0
    t1 = 200.0
    t2 = 200.0
    t3 = 0.0
    t4 = 0.0

    u0 = initializeField(imax, jmax, t0, t1, t2, t3, t4)

    if do_sim or do_all:
        sim(
            u0, t_end, deltax, deltay,
            alpha, dt_implicit, dt_explicit_given,
            t1, t2, t3, t4
        )

    if do_conv or do_all:
        convergence(
            u0, t_end, deltax, deltay,
            alpha, t1, t2, t3, t4
        )

    print("\nDone")


if __name__ == "__main__":
    main()