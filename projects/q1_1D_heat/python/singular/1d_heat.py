# ============================================================================
# Heat equation 1D (multiple schemes)  Python translation of the C++ layout
# Uses 1-based indexing convention: valid indices are 1..imax, index 0 unused
# ============================================================================
import os
import math
import numpy as np
import matplotlib.pyplot as plt


# ============================================================================
# Helpers for directories
# ============================================================================
def ensure_dir(path: str) -> None:
    os.makedirs(path, exist_ok=True)


def ensure_plot_dir() -> None:
    ensure_dir("plot")


def ensure_compare_dir() -> None:
    ensure_dir("compare")


def ensure_data_dir() -> None:
    ensure_dir("data")


def make_dt_tag(dt: float) -> str:
    code = int(round(dt * 100.0))
    return f"{code:03d}"


def file_exists(path: str) -> bool:
    return os.path.isfile(path)


def scheme_data_exists(scheme: str, dt: float) -> bool:
    tag = make_dt_tag(dt)
    fname = os.path.join("data", f"{scheme}_{tag}.txt")
    return file_exists(fname)


def last_block_index_and_time(path: str):
    if not file_exists(path):
        return 0, 0.0

    idx = -1
    last_t = 0.0
    with open(path, "r", encoding="utf-8") as f:
        for line in f:
            if line.startswith("# t ="):
                idx += 1
                try:
                    last_t = float(line.split("=")[1].split()[0])
                except Exception:
                    pass

    if idx < 0:
        idx = 0
    return idx, last_t


def latest_time_in_file(path: str) -> float:
    if not file_exists(path):
        return 0.0
    last_t = 0.0
    with open(path, "r", encoding="utf-8") as f:
        for line in f:
            if line.startswith("# t ="):
                try:
                    last_t = float(line.split("=")[1].split()[0])
                except Exception:
                    pass
    return last_t


# ============================================================================
# Error norms
# ============================================================================
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


# ============================================================================
# Plotting and comparison (matplotlib versions instead of gnuplot)
# ============================================================================
def load_block_from_file(path: str, idx: int):
    if not file_exists(path):
        return None, None

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
                    pass
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
        return None, None
    return current_t, (np.array(x), np.array(y))

def plot_one_dt(delta_t: float) -> None:
    ensure_plot_dir()
    tag = make_dt_tag(delta_t)

    # --------------------------------------------------
    # Files
    # --------------------------------------------------
    exact_file = os.path.join("data", f"exact_{tag}.txt")
    cn_file    = os.path.join("data", f"cn_{tag}.txt")
    ftcs_e_file = os.path.join("data", f"ftcs_explicit_{tag}.txt")
    dufort_file = os.path.join("data", f"dufort_{tag}.txt")
    ftcs_i_file = os.path.join("data", f"ftcs_implicit_{tag}.txt")

    # --------------------------------------------------
    # Time index
    # --------------------------------------------------
    idx, tlabel = last_block_index_and_time(exact_file)
    if tlabel == 0.0:
        idx = 4
        tlabel = 0.4

    # --------------------------------------------------
    # Load data
    # --------------------------------------------------
    _, data = load_block_from_file(exact_file, idx)
    x, exact = data if data is not None else (None, None)

    _, data = load_block_from_file(cn_file, idx)
    x, cn = data if data is not None else (None, None)

    _, data = load_block_from_file(ftcs_e_file, idx)
    x, ftcs_explicit = data if data is not None else (None, None)

    _, data = load_block_from_file(dufort_file, idx)
    x, dufort = data if data is not None else (None, None)

    _, data = load_block_from_file(ftcs_i_file, idx)
    x, ftcs_implicit = data if data is not None else (None, None)

    # --------------------------------------------------
    # Plot
    # --------------------------------------------------
    plt.figure()
    plt.title(f"All schemes, dt={delta_t:.3f}, t={tlabel:.2f}")
    plt.xlabel("X")
    plt.ylabel("T")
    plt.grid(True)

    if exact is not None:
        plt.plot(
            x, exact,
            color="black",
            linestyle="--",
            linewidth=3.0,
            label="exact"
        )

    if cn is not None:
        plt.plot(
            x, cn,
            marker="d",
            linestyle="-",
            linewidth=1.5,
            label="cn"
        )

    if ftcs_explicit is not None:
        plt.plot(
            x, ftcs_explicit,
            marker="o",
            linestyle="-",
            linewidth=1.5,
            label="ftcs explicit"
        )

    if dufort is not None:
        plt.plot(
            x, dufort,
            marker="s",
            linestyle="-",
            linewidth=1.5,
            label="dufort"
        )

    if ftcs_implicit is not None:
        plt.plot(
            x, ftcs_implicit,
            marker="^",
            linestyle="-",
            linewidth=1.5,
            label="ftcs implicit"
        )

    plt.legend()

    outpng = os.path.join("plot", f"all_schemes_{tag}.png")
    plt.savefig(outpng, dpi=150)
    plt.show()



def compare_error_schemes(delta_t: float, t_target: float = None) -> None:
    ensure_compare_dir()
    tag = make_dt_tag(delta_t)

    exact_path = os.path.join("data", f"exact_{tag}.txt")
    if t_target is None:
        t_latest = latest_time_in_file(exact_path)
        if t_latest <= 0.0:
            return
        t_target = t_latest

    idx = int(round(t_target / 0.1))

    schemes = [
        ("ftcs_explicit", os.path.join("data", f"ftcs_explicit_{tag}.txt")),
        ("dufort",        os.path.join("data", f"dufort_{tag}.txt")),
        ("ftcs_implicit", os.path.join("data", f"ftcs_implicit_{tag}.txt")),
        ("cn",            os.path.join("data", f"cn_{tag}.txt")),
    ]

    t_ex, data_ex = load_block_from_file(exact_path, idx)
    if data_ex is None:
        return
    x_ex, y_ex = data_ex

    plt.figure()
    for name, path in schemes:
        t_s, data_s = load_block_from_file(path, idx)
        if data_s is None:
            continue
        x_s, y_s = data_s
        if len(x_s) != len(x_ex):
            continue
        plt.plot(x_s, y_s - y_ex, marker="o", label=name)

    plt.title(f"Errors vs exact, dt={delta_t:.3f}, t={t_target:.2f}")
    plt.xlabel("x")
    plt.ylabel("T scheme minus T exact")
    plt.grid(True)
    plt.legend()
    outpng = os.path.join("compare", f"error_schemes_{tag}_t{t_target:.1f}.png")
    plt.savefig(outpng, dpi=150)


def convergence_study(delta_x: float, imax: int, alpha: float, t0: float,
                      tboundary: float, t_target: float,
                      dt_list: list[float]) -> None:
    ensure_data_dir()
    ensure_compare_dir()

    outdata = os.path.join("data", f"convergence_t{int(round(t_target * 100.0)):03d}.txt")

    with open(outdata, "w", encoding="utf-8") as f:
        f.write(f"# convergence at t = {t_target:.2f} hr\n")
        f.write("# dt\tL2_exp\tL2_imp\tL2_df\tL2_cn\tLinf_exp\tLinf_imp\tLinf_df\tLinf_cn\n")

        x_vector = np.zeros(imax + 1, dtype=float)
        x_vector[1] = 0.0
        for i in range(1, imax):
            x_vector[i + 1] = x_vector[i] + delta_x

        L = (imax - 1) * delta_x

        for dt in dt_list:
            F_num = (alpha * dt) / (delta_x * delta_x)
            nmax = int(round(t_target / dt))

            T_exp = np.zeros(imax + 1, dtype=float)
            T_imp = np.zeros(imax + 1, dtype=float)
            T_df = np.zeros(imax + 1, dtype=float)
            T_cn = np.zeros(imax + 1, dtype=float)
            T_ex = np.zeros(imax + 1, dtype=float)

            FTCS_explicit(nmax, F_num, tboundary, t0, imax, T_exp)
            FTCS_implicit(nmax, F_num, tboundary, t0, imax, T_imp)
            DufortFrankel(nmax, F_num, tboundary, t0, imax, T_df)
            CrankNicolson(nmax, F_num, tboundary, t0, imax, T_cn)

            exact_solution_profile(imax, x_vector, t_target, alpha, L, tboundary, t0, 200, T_ex)

            L2_exp = err_L2(imax, T_exp, T_ex)
            L2_imp = err_L2(imax, T_imp, T_ex)
            L2_df = err_L2(imax, T_df, T_ex)
            L2_cn = err_L2(imax, T_cn, T_ex)

            Li_exp = err_Linf(imax, T_exp, T_ex)
            Li_imp = err_Linf(imax, T_imp, T_ex)
            Li_df = err_Linf(imax, T_df, T_ex)
            Li_cn = err_Linf(imax, T_cn, T_ex)

            f.write(
                f"{dt:.6f}\t{L2_exp:.10e}\t{L2_imp:.10e}\t{L2_df:.10e}\t{L2_cn:.10e}\t"
                f"{Li_exp:.10e}\t{Li_imp:.10e}\t{Li_df:.10e}\t{Li_cn:.10e}\n"
            )

    data = np.loadtxt(outdata, comments="#")
    dt = data[:, 0]
    L2_exp, L2_imp, L2_df, L2_cn = data[:, 1], data[:, 2], data[:, 3], data[:, 4]
    Li_exp, Li_imp, Li_df, Li_cn = data[:, 5], data[:, 6], data[:, 7], data[:, 8]

    plt.figure()
    plt.loglog(dt, L2_exp, marker="o", label="L2 exp")
    plt.loglog(dt, L2_imp, marker="o", label="L2 imp")
    plt.loglog(dt, L2_df, marker="o", label="L2 df")
    plt.loglog(dt, L2_cn, marker="o", label="L2 cn")
    plt.gca().invert_xaxis()
    plt.title(f"Convergence at t={t_target:.2f}")
    plt.xlabel("dt")
    plt.ylabel("L2 error")
    plt.grid(True, which="both")
    plt.legend()
    outpng = os.path.join("compare", f"convergence_t{t_target:.2f}.png")
    plt.savefig(outpng, dpi=150)


# ============================================================================
# Simulation driver (writes data files like the C++ sim)
# ============================================================================
def sim(delta_x: float, delta_t: float, imax: int,
        t0: float, tboundary: float, F_num: float) -> None:
    ensure_data_dir()

    tag = make_dt_tag(delta_t)

    f_explicit = os.path.join("data", f"ftcs_explicit_{tag}.txt")
    f_dufort = os.path.join("data", f"dufort_{tag}.txt")
    f_implicit = os.path.join("data", f"ftcs_implicit_{tag}.txt")
    f_cn = os.path.join("data", f"cn_{tag}.txt")
    f_exact = os.path.join("data", f"exact_{tag}.txt")

    f_err_exp = os.path.join("data", f"error_ftcs_explicit_{tag}.txt")
    f_err_imp = os.path.join("data", f"error_ftcs_implicit_{tag}.txt")
    f_err_df = os.path.join("data", f"error_dufort_{tag}.txt")
    f_err_cn = os.path.join("data", f"error_cn_{tag}.txt")

    x_vector = np.zeros(imax + 1, dtype=float)
    x_vector[1] = 0.0
    for i in range(1, imax):
        x_vector[i + 1] = x_vector[i] + delta_x

    L = (imax - 1) * delta_x
    alpha = F_num * (delta_x * delta_x) / delta_t

    with open(f_explicit, "w", encoding="utf-8") as out_exp, \
         open(f_dufort, "w", encoding="utf-8") as out_df, \
         open(f_implicit, "w", encoding="utf-8") as out_imp, \
         open(f_cn, "w", encoding="utf-8") as out_cn, \
         open(f_exact, "w", encoding="utf-8") as out_ex, \
         open(f_err_exp, "w", encoding="utf-8") as out_eexp, \
         open(f_err_imp, "w", encoding="utf-8") as out_eimp, \
         open(f_err_df, "w", encoding="utf-8") as out_edf, \
         open(f_err_cn, "w", encoding="utf-8") as out_ecn:

        for k in range(0, 5):
            t_target = 0.1 * k
            nmax = int(round(t_target / delta_t))

            T_exact = np.zeros(imax + 1, dtype=float)
            T_explicit = np.zeros(imax + 1, dtype=float)
            T_dufort = np.zeros(imax + 1, dtype=float)
            T_implicit = np.zeros(imax + 1, dtype=float)
            T_cn = np.zeros(imax + 1, dtype=float)

            FTCS_explicit(nmax, F_num, tboundary, t0, imax, T_explicit)
            DufortFrankel(nmax, F_num, tboundary, t0, imax, T_dufort)
            FTCS_implicit(nmax, F_num, tboundary, t0, imax, T_implicit)
            CrankNicolson(nmax, F_num, tboundary, t0, imax, T_cn)

            exact_solution_profile(imax, x_vector, t_target, alpha, L, tboundary, t0, 200, T_exact)

            header = f"# t = {t_target:.2f} hr\n# x\tT\n"
            out_exp.write(header)
            out_df.write(header)
            out_imp.write(header)
            out_cn.write(header)
            out_ex.write(header)

            eheader = f"# t = {t_target:.2f} hr\n# x\terr\n"
            out_eexp.write(eheader)
            out_eimp.write(eheader)
            out_edf.write(eheader)
            out_ecn.write(eheader)

            for i in range(1, imax + 1):
                out_exp.write(f"{x_vector[i]:.6f}\t{T_explicit[i]:.6f}\n")
                out_df.write(f"{x_vector[i]:.6f}\t{T_dufort[i]:.6f}\n")
                out_imp.write(f"{x_vector[i]:.6f}\t{T_implicit[i]:.6f}\n")
                out_cn.write(f"{x_vector[i]:.6f}\t{T_cn[i]:.6f}\n")
                out_ex.write(f"{x_vector[i]:.6f}\t{T_exact[i]:.6f}\n")

                out_eexp.write(f"{x_vector[i]:.6f}\t{(T_explicit[i] - T_exact[i]):.6f}\n")
                out_eimp.write(f"{x_vector[i]:.6f}\t{(T_implicit[i] - T_exact[i]):.6f}\n")
                out_edf.write(f"{x_vector[i]:.6f}\t{(T_dufort[i] - T_exact[i]):.6f}\n")
                out_ecn.write(f"{x_vector[i]:.6f}\t{(T_cn[i] - T_exact[i]):.6f}\n")

            out_exp.write("\n\n")
            out_df.write("\n\n")
            out_imp.write("\n\n")
            out_cn.write("\n\n")
            out_ex.write("\n\n")

            out_eexp.write("\n\n")
            out_eimp.write("\n\n")
            out_edf.write("\n\n")
            out_ecn.write("\n\n")


# ============================================================================
# Main, mirrors the C++ main
# ============================================================================
def main():
    imax = 21
    delta_x = 0.05
    delta_t = 0.01
    alpha = 0.1
    t0 = 100.0
    tboundary = 300.0

    Fourier_num = (alpha * delta_t) / (delta_x ** 2.0)

    sim(delta_x, delta_t, imax, t0, tboundary, Fourier_num)

    print("Plotting all schemes at 1 dt\n")
    plot_one_dt(delta_t)

    # print("Comparing different schemes\n")
    # compare_error_schemes(delta_t, 0.4)

    # dt_list = [0.10, 0.05, 0.02, 0.01, 0.005]
    # print("Convergence study\n")
    # convergence_study(delta_x, imax, alpha, t0, tboundary, 0.4, dt_list)


if __name__ == "__main__":
    main()
