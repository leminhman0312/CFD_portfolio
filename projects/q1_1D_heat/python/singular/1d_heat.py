# ============================================================================
# heat_1d_clean.py
#
# Learning-focused 1D heat equation code
#
# Structure
#   - main
#   - helpers
#   - numerical schemes
#   - sim (writes data files)
#   - convergence study (writes one table)
#   - plot (reads files and makes pngs)
#
# Run
#   python3 heat_1d_clean.py
# ============================================================================

import os
import math
import numpy as np
import matplotlib.pyplot as plt


# ============================================================================
# Helpers
# ============================================================================

def create_dir() -> None:
    os.makedirs("data", exist_ok=True)
    os.makedirs("plot", exist_ok=True)


def make_dt_tag(dt: float) -> str:
    code = int(round(dt * 100.0))
    return f"{code:03d}"


def make_t_tag(t: float) -> str:
    code = int(round(t * 100.0))
    return f"{code:03d}"


def build_grid(imax: int, dx: float) -> np.ndarray:
    x = np.zeros(imax + 1, dtype=float)
    x[1] = 0.0
    for i in range(1, imax):
        x[i + 1] = x[i] + dx
    return x


def load_block_from_file(path: str, idx: int):
    if not os.path.isfile(path):
        return None, None, False

    block = -1
    x = []
    y = []
    current_t = 0.0

    with open(path, "r", encoding="utf-8") as f:
        for line in f:
            s = line.strip()
            if not s:
                continue

            if s.startswith("# t ="):
                block += 1
                if block > idx:
                    break
                try:
                    current_t = float(s.split("=")[1].split()[0])
                except Exception:
                    current_t = 0.0
                x = []
                y = []
                continue

            if s.startswith("#"):
                continue

            if block == idx:
                parts = s.split()
                if len(parts) >= 2:
                    x.append(float(parts[0]))
                    y.append(float(parts[1]))

    if block < idx:
        return None, None, False

    return current_t, (np.array(x), np.array(y)), True


def write_block_T(f, t: float, x: np.ndarray, T: np.ndarray, imax: int) -> None:
    f.write(f"# t = {t:.2f} hr\n")
    f.write("# x\tT\n")
    for i in range(1, imax + 1):
        f.write(f"{x[i]:.6f}\t{T[i]:.6f}\n")
    f.write("\n\n")


def write_block_error(f, t: float, x: np.ndarray, Terr: np.ndarray, imax: int) -> None:
    f.write(f"# t = {t:.2f} hr\n")
    f.write("# x\terr\n")
    for i in range(1, imax + 1):
        f.write(f"{x[i]:.6f}\t{Terr[i]:.6f}\n")
    f.write("\n\n")


def err_L2(imax: int, a: np.ndarray, b: np.ndarray) -> float:
    s = 0.0
    for i in range(1, imax + 1):
        e = a[i] - b[i]
        s += e * e
    return math.sqrt(s / float(imax))


def err_Linf(imax: int, a: np.ndarray, b: np.ndarray) -> float:
    m = 0.0
    for i in range(1, imax + 1):
        e = abs(a[i] - b[i])
        if e > m:
            m = e
    return m


# ============================================================================
# Thomas tridiagonal solver
# Solves u[1..N] for arrays a,b,c,d indexed 1..N
# ============================================================================

def thomasTriDiagonal(N: int, a: np.ndarray, b: np.ndarray, c: np.ndarray,
                      d: np.ndarray, u: np.ndarray) -> None:
    dprime = np.zeros(N + 1, dtype=float)
    cprime = np.zeros(N + 1, dtype=float)

    dprime[1] = d[1]
    cprime[1] = c[1]

    for i in range(2, N + 1):
        dprime[i] = d[i] - (b[i] * a[i - 1]) / dprime[i - 1]
        cprime[i] = c[i] - (cprime[i - 1] * b[i]) / dprime[i - 1]

    u[N] = cprime[N] / dprime[N]
    for i in range(N - 1, 0, -1):
        u[i] = (cprime[i] - a[i] * u[i + 1]) / dprime[i]


# ============================================================================
# Exact solution (series)
# ============================================================================

def exact_solution_profile(imax: int, x: np.ndarray, t: float,
                           alpha: float, L: float,
                           tb: float, t0: float,
                           nterms: int, Tout: np.ndarray) -> None:
    pi = 4.0 * math.atan(1.0)
    A = t0 - tb

    for i in range(1, imax + 1):
        s = 0.0
        xi = x[i]
        for n in range(1, nterms + 1):
            coeff = (2.0 * A / (n * pi)) * (1.0 - ((-1.0) ** n))
            if coeff == 0.0:
                continue
            k = n * pi / L
            s += coeff * math.sin(k * xi) * math.exp(-alpha * k * k * t)
        Tout[i] = tb + s

    Tout[1] = tb
    Tout[imax] = tb


# ============================================================================
# Numerical schemes
# All arrays length imax+1 and valid indices are 1..imax
# ============================================================================

def FTCS_explicit(nmax: int, F: float, tb: float, t0: float,
                  imax: int, Tout: np.ndarray) -> None:
    u = np.zeros(imax + 1, dtype=float)
    un = np.zeros(imax + 1, dtype=float)

    un[1] = tb
    un[imax] = tb
    for i in range(2, imax):
        un[i] = t0

    for _ in range(1, nmax + 1):
        u[:] = un[:]
        for i in range(2, imax):
            un[i] = u[i] + F * (u[i + 1] - 2.0 * u[i] + u[i - 1])
        un[1] = tb
        un[imax] = tb

    Tout[:] = un[:]


def implicit_interior_one_step(F: float, tb: float, t0: float,
                               imax: int, u_full: np.ndarray) -> None:
    N = imax - 2
    a = np.zeros(N + 1, dtype=float)
    b = np.zeros(N + 1, dtype=float)
    c = np.zeros(N + 1, dtype=float)
    d = np.zeros(N + 1, dtype=float)
    u_int = np.zeros(N + 1, dtype=float)

    u_full[1] = tb
    u_full[imax] = tb
    for j in range(1, N + 1):
        u_full[j + 1] = t0

    for j in range(1, N + 1):
        d[j] = 1.0 + 2.0 * F
        a[j] = -F
        b[j] = -F

    b[1] = 0.0
    a[N] = 0.0

    for j in range(1, N + 1):
        c[j] = u_full[j + 1]

    c[1] += F * tb
    c[N] += F * tb

    thomasTriDiagonal(N, a, b, c, d, u_int)

    u_full[1] = tb
    u_full[imax] = tb
    for j in range(1, N + 1):
        u_full[j + 1] = u_int[j]


def DufortFrankel(nmax: int, F: float, tb: float, t0: float,
                  imax: int, Tout: np.ndarray) -> None:
    dval = 2.0 * F

    un_m1 = np.zeros(imax + 1, dtype=float)
    un = np.zeros(imax + 1, dtype=float)
    un_p1 = np.zeros(imax + 1, dtype=float)
    u1 = np.zeros(imax + 1, dtype=float)

    un_m1[1] = tb
    un_m1[imax] = tb
    for i in range(2, imax):
        un_m1[i] = t0

    implicit_interior_one_step(F, tb, t0, imax, u1)
    un[:] = u1[:]

    for _ in range(1, nmax):
        for i in range(2, imax):
            un_p1[i] = ((1.0 - dval) * un_m1[i] + dval * (un[i + 1] + un[i - 1])) / (1.0 + dval)

        un_p1[1] = tb
        un_p1[imax] = tb

        un_m1[:] = un[:]
        un[:] = un_p1[:]

    Tout[:] = un[:]


def FTCS_implicit(nmax: int, F: float, tb: float, t0: float,
                  imax: int, Tout: np.ndarray) -> None:
    a = np.zeros(imax + 1, dtype=float)
    b = np.zeros(imax + 1, dtype=float)
    c = np.zeros(imax + 1, dtype=float)
    d = np.zeros(imax + 1, dtype=float)
    un = np.zeros(imax + 1, dtype=float)

    un[1] = tb
    un[imax] = tb
    for i in range(2, imax):
        un[i] = t0

    d[1] = 1.0
    d[imax] = 1.0

    for i in range(2, imax):
        d[i] = 1.0 + 2.0 * F
        a[i] = -F
        b[i] = -F

    for _ in range(1, nmax + 1):
        c[1] = tb
        c[imax] = tb
        for i in range(2, imax):
            c[i] = un[i]

        thomasTriDiagonal(imax, a, b, c, d, un)

        un[1] = tb
        un[imax] = tb

    Tout[:] = un[:]


def CrankNicolson(nmax: int, F: float, tb: float, t0: float,
                  imax: int, Tout: np.ndarray) -> None:
    dnc = F / 2.0

    a = np.zeros(imax + 1, dtype=float)
    b = np.zeros(imax + 1, dtype=float)
    c = np.zeros(imax + 1, dtype=float)
    d = np.zeros(imax + 1, dtype=float)
    u0 = np.zeros(imax + 1, dtype=float)
    uhalf = np.zeros(imax + 1, dtype=float)

    u0[1] = tb
    u0[imax] = tb
    for i in range(2, imax):
        u0[i] = t0

    d[1] = 1.0
    d[imax] = 1.0

    for i in range(2, imax):
        d[i] = 1.0 + 2.0 * dnc
        a[i] = -dnc
        b[i] = -dnc

    for _ in range(1, nmax + 1):
        uhalf[1] = tb
        uhalf[imax] = tb
        for i in range(2, imax):
            uhalf[i] = u0[i] + dnc * (u0[i + 1] - 2.0 * u0[i] + u0[i - 1])

        c[1] = tb
        c[imax] = tb
        for i in range(2, imax):
            c[i] = uhalf[i]

        thomasTriDiagonal(imax, a, b, c, d, u0)

        u0[1] = tb
        u0[imax] = tb

    Tout[:] = u0[:]


# ============================================================================
# sim: write profile and error files at t = 0.0 0.1 0.2 0.3 0.4
# ============================================================================

def sim(dx: float, dt: float, imax: int, t0: float, tb: float, F: float) -> None:
    tag = make_dt_tag(dt)

    fexp = os.path.join("data", f"ftcs_explicit_{tag}.txt")
    fdf  = os.path.join("data", f"dufort_{tag}.txt")
    fimp = os.path.join("data", f"ftcs_implicit_{tag}.txt")
    fcn  = os.path.join("data", f"cn_{tag}.txt")
    fex  = os.path.join("data", f"exact_{tag}.txt")

    ferr_exp = os.path.join("data", f"error_ftcs_explicit_{tag}.txt")
    ferr_df  = os.path.join("data", f"error_dufort_{tag}.txt")
    ferr_imp = os.path.join("data", f"error_ftcs_implicit_{tag}.txt")
    ferr_cn  = os.path.join("data", f"error_cn_{tag}.txt")

    x = build_grid(imax, dx)
    L = (imax - 1) * dx
    alpha = F * (dx * dx) / dt

    Tex = np.zeros(imax + 1, dtype=float)
    Texp = np.zeros(imax + 1, dtype=float)
    Tdf = np.zeros(imax + 1, dtype=float)
    Timp = np.zeros(imax + 1, dtype=float)
    Tcn = np.zeros(imax + 1, dtype=float)

    with open(fexp, "w", encoding="utf-8") as out_exp, \
         open(fdf,  "w", encoding="utf-8") as out_df, \
         open(fimp, "w", encoding="utf-8") as out_imp, \
         open(fcn,  "w", encoding="utf-8") as out_cn, \
         open(fex,  "w", encoding="utf-8") as out_ex, \
         open(ferr_exp, "w", encoding="utf-8") as out_eexp, \
         open(ferr_df,  "w", encoding="utf-8") as out_edf, \
         open(ferr_imp, "w", encoding="utf-8") as out_eimp, \
         open(ferr_cn,  "w", encoding="utf-8") as out_ecn:

        for k in range(0, 5):
            t = 0.1 * k
            nmax = int(round(t / dt))

            FTCS_explicit(nmax, F, tb, t0, imax, Texp)
            DufortFrankel(nmax, F, tb, t0, imax, Tdf)
            FTCS_implicit(nmax, F, tb, t0, imax, Timp)
            CrankNicolson(nmax, F, tb, t0, imax, Tcn)

            exact_solution_profile(imax, x, t, alpha, L, tb, t0, 200, Tex)

            write_block_T(out_exp, t, x, Texp, imax)
            write_block_T(out_df,  t, x, Tdf,  imax)
            write_block_T(out_imp, t, x, Timp, imax)
            write_block_T(out_cn,  t, x, Tcn,  imax)
            write_block_T(out_ex,  t, x, Tex,  imax)

            write_block_error(out_eexp, t, x, Texp - Tex, imax)
            write_block_error(out_edf,  t, x, Tdf  - Tex, imax)
            write_block_error(out_eimp, t, x, Timp - Tex, imax)
            write_block_error(out_ecn,  t, x, Tcn  - Tex, imax)


# ============================================================================
# convergence: write one convergence table at t_target
# ============================================================================

def convergence_study(dx: float, imax: int, alpha: float,
                      t0: float, tb: float, t_target: float) -> str:
    dt_list = [0.10, 0.05, 0.02, 0.01, 0.005]
    outdata = os.path.join("data", f"convergence_t{make_t_tag(t_target)}.txt")

    x = build_grid(imax, dx)
    L = (imax - 1) * dx

    Tex = np.zeros(imax + 1, dtype=float)
    Texp = np.zeros(imax + 1, dtype=float)
    Timp = np.zeros(imax + 1, dtype=float)
    Tdf = np.zeros(imax + 1, dtype=float)
    Tcn = np.zeros(imax + 1, dtype=float)

    with open(outdata, "w", encoding="utf-8") as f:
        f.write(f"# convergence at t = {t_target:.2f} hr\n")
        f.write("# dt\tL2_exp\tL2_imp\tL2_df\tL2_cn\tLinf_exp\tLinf_imp\tLinf_df\tLinf_cn\n")

        for dt in dt_list:
            F = alpha * dt / (dx * dx)
            nmax = int(round(t_target / dt))

            FTCS_explicit(nmax, F, tb, t0, imax, Texp)
            FTCS_implicit(nmax, F, tb, t0, imax, Timp)
            DufortFrankel(nmax, F, tb, t0, imax, Tdf)
            CrankNicolson(nmax, F, tb, t0, imax, Tcn)

            exact_solution_profile(imax, x, t_target, alpha, L, tb, t0, 200, Tex)

            L2_exp = err_L2(imax, Texp, Tex)
            L2_imp = err_L2(imax, Timp, Tex)
            L2_df  = err_L2(imax, Tdf,  Tex)
            L2_cn  = err_L2(imax, Tcn,  Tex)

            Li_exp = err_Linf(imax, Texp, Tex)
            Li_imp = err_Linf(imax, Timp, Tex)
            Li_df  = err_Linf(imax, Tdf,  Tex)
            Li_cn  = err_Linf(imax, Tcn,  Tex)

            f.write(
                f"{dt:.6f}\t{L2_exp:.10e}\t{L2_imp:.10e}\t{L2_df:.10e}\t{L2_cn:.10e}\t"
                f"{Li_exp:.10e}\t{Li_imp:.10e}\t{Li_df:.10e}\t{Li_cn:.10e}\n"
            )

    return outdata


# ============================================================================
# plot: make all pngs for one dt plus convergence
# ============================================================================

def plot(dt: float, t_target: float) -> None:
    tag = make_dt_tag(dt)

    # blocks in files are:
    # idx 0 -> t=0.0, idx 1 -> t=0.1, idx 2 -> t=0.2, idx 3 -> t=0.3, idx 4 -> t=0.4
    idx_target = int(round(t_target / 0.1))
    if idx_target < 0:
        idx_target = 0

    # -------------------------------
    # helper: plot one scheme file (4 time blocks)
    # -------------------------------
    def plot_scheme_times(scheme: str, title_name: str, outpng: str) -> None:
        path = os.path.join("data", f"{scheme}_{tag}.txt")

        plt.figure()
        plt.xlabel("X [ft]")
        plt.ylabel("T [degree F]")
        plt.title(f"∂T/∂t = α ∂²T/∂x² , {title_name}, dt = {dt:.2f} hr")
        plt.grid(True)
        plt.xlim(0.0, 1.0)
        plt.ylim(100.0, 300.0)

        # we want blocks 1..4 (t=0.1..0.4), not block 0
        for block_idx, tlabel in [(1, "0.1"), (2, "0.2"), (3, "0.3"), (4, "0.4")]:
            t_block, data, ok = load_block_from_file(path, block_idx)
            if not ok or data is None:
                continue
            x, y = data
            plt.plot(x, y, marker="o", linewidth=2.0, label=f"t = {tlabel} hr")

        plt.legend()
        plt.savefig(outpng, dpi=150, bbox_inches="tight")
        plt.close()

    # -------------------------------
    # 1..4 scheme plots
    # -------------------------------
    plot_scheme_times("ftcs_explicit", "FTCS explicit", os.path.join("plot", f"ftcs_explicit_{tag}.png"))
    plot_scheme_times("dufort",        "Dufort Frankel", os.path.join("plot", f"dufort_{tag}.png"))
    plot_scheme_times("ftcs_implicit", "FTCS implicit", os.path.join("plot", f"ftcs_implicit_{tag}.png"))
    plot_scheme_times("cn",            "Crank Nicolson", os.path.join("plot", f"cn_{tag}.png"))

    # -------------------------------
    # 5) errors vs exact at t_target
    # -------------------------------
    exact_file = os.path.join("data", f"exact_{tag}.txt")
    t_ex, data_ex, ok_ex = load_block_from_file(exact_file, idx_target)
    if not ok_ex or data_ex is None:
        raise RuntimeError("Could not load exact block at t_target. Run sim first.")

    x_ex, y_ex = data_ex

    plt.figure()
    plt.xlabel("X [ft]")
    plt.ylabel("Error = T_scheme - T_exact [°F]")
    plt.title(f"Error vs Exact, dt = {dt:.2f} hr, t = {t_target:.1f} hr")
    plt.grid(True)
    plt.xlim(0.0, 1.0)

    # exact error line at 0
    plt.plot(x_ex, np.zeros_like(x_ex), linewidth=2.5, linestyle="--", label="Exact error (0)")

    for scheme, name in [
        ("ftcs_explicit", "FTCS explicit"),
        ("ftcs_implicit", "FTCS implicit"),
        ("dufort",        "Dufort Frankel"),
        ("cn",            "Crank Nicolson"),
    ]:
        path = os.path.join("data", f"{scheme}_{tag}.txt")
        _, data_s, ok_s = load_block_from_file(path, idx_target)
        if not ok_s or data_s is None:
            continue
        xs, ys = data_s

        # if x grids differ slightly, interpolate scheme onto exact x grid
        if len(xs) != len(x_ex) or np.max(np.abs(xs - x_ex)) > 1.0e-12:
            ys_i = np.interp(x_ex, xs, ys)
            err = ys_i - y_ex
            plt.plot(x_ex, err, marker="o", linewidth=1.5, label=name)
        else:
            err = ys - y_ex
            plt.plot(xs, err, marker="o", linewidth=1.5, label=name)

    plt.legend()
    out_err = os.path.join("plot", f"error_schemes_{tag}_t{t_target:.1f}.png")
    plt.savefig(out_err, dpi=150, bbox_inches="tight")
    plt.close()

    # -------------------------------
    # 6) convergence plot from table (L2 only, like your gnuplot)
    # -------------------------------
    conv_path = os.path.join("data", f"convergence_t{make_t_tag(t_target)}.txt")
    if os.path.isfile(conv_path):
        data = np.loadtxt(conv_path, comments="#")
        dtc = data[:, 0]
        L2_exp, L2_imp, L2_df, L2_cn = data[:, 1], data[:, 2], data[:, 3], data[:, 4]

        plt.figure()
        plt.xlabel("dt [hr]")
        plt.ylabel("Error norm")
        plt.title(f"Convergence vs dt at t = {t_target:.2f} hr")
        plt.grid(True, which="both")
        plt.xscale("log")
        plt.yscale("log")

        plt.plot(dtc, L2_exp, marker="o", linewidth=2.0, label="FTCS explicit L2")
        plt.plot(dtc, L2_imp, marker="o", linewidth=2.0, label="FTCS implicit L2")
        plt.plot(dtc, L2_df,  marker="o", linewidth=2.0, label="Dufort Frankel L2")
        plt.plot(dtc, L2_cn,  marker="o", linewidth=2.0, label="Crank Nicolson L2")

        plt.legend()
        out_conv = os.path.join("plot", f"convergence_t{t_target:.2f}.png")
        plt.savefig(out_conv, dpi=150, bbox_inches="tight")
        plt.close()


# ============================================================================
# main
# ============================================================================

def main():
    imax = 21
    dx = 0.05
    dt = 0.01
    alpha = 0.1
    t0 = 100.0
    tb = 300.0

    t_target = 0.4
    F = alpha * dt / (dx * dx)

    create_dir()
    sim(dx, dt, imax, t0, tb, F)
    convergence_study(dx, imax, alpha, t0, tb, t_target)
    plot(dt, t_target)

    print("Done")


if __name__ == "__main__":
    main()
