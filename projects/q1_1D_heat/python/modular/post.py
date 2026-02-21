import os
import numpy as np
import matplotlib.pyplot as plt

from io_utils import (
    create_dir,
    make_dt_tag,
    make_t_tag,
    load_block_from_file,
    write_block_T,
    write_block_error
)

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
