#!/usr/bin/env python3

import os
import sys
import numpy as np
import matplotlib.pyplot as plt


def make_dt_tag(dt):
    code = round(dt * 100.0)
    return f"{code:03d}"


def make_t_tag(t):
    code = round(t * 100.0)
    return f"{code:03d}"


def scheme_title(scheme):
    if scheme == "ftcs_explicit":
        return "FTCS explicit"
    if scheme == "ftcs_implicit":
        return "FTCS implicit"
    if scheme == "dufort":
        return "Dufort-Frankel"
    if scheme == "cn":
        return "Crank-Nicolson"
    return scheme


def read_blocks(filename):
    """
    Reads files of the form:

    # t = 0.10 hr
    # x   T
    ...
    <blank lines>

    Returns a list of dicts:
    [
      {"t": 0.0, "x": array(...), "y": array(...)},
      ...
    ]
    """
    blocks = []
    current_t = None
    xvals = []
    yvals = []

    def flush_block():
        nonlocal current_t, xvals, yvals
        if current_t is not None and len(xvals) > 0:
            blocks.append({
                "t": current_t,
                "x": np.array(xvals, dtype=float),
                "y": np.array(yvals, dtype=float),
            })
        current_t = None
        xvals = []
        yvals = []

    with open(filename, "r") as f:
        for raw in f:
            line = raw.strip()

            if not line:
                flush_block()
                continue

            if line.startswith("# t ="):
                parts = line.split()
                # "# t = 0.10 hr"
                current_t = float(parts[3])
                continue

            if line.startswith("#"):
                continue

            parts = line.split()
            if len(parts) >= 2:
                xvals.append(float(parts[0]))
                yvals.append(float(parts[1]))

    flush_block()
    return blocks


def ensure_plot_dir():
    os.makedirs("plot", exist_ok=True)


def plot_scheme(dt, scheme, tag):
    """
    Mirrors plot_scheme.gp

    Original gnuplot used:
      index 1, 2, 3, 4
    which correspond to t = 0.1, 0.2, 0.3, 0.4
    """
    fname = f"data/{scheme}_{tag}.txt"
    blocks = read_blocks(fname)

    fig, ax = plt.subplots(figsize=(10, 7))

    markers = ["o", "s", "^", "d"]
    target_indices = [1, 2, 3, 4]

    for k, idx in enumerate(target_indices):
        if idx < len(blocks):
            b = blocks[idx]
            ax.plot(
                b["x"], b["y"],
                marker=markers[k],
                linewidth=2,
                markersize=6,
                label=f"t = {b['t']:.1f} hr"
            )

    ax.set_xlabel("X [ft]")
    ax.set_ylabel("T [degree F]")
    ax.set_title(
        f"∂T/∂t = α [∂²T/∂x²] , {scheme_title(scheme)}, dt = {dt:.2f} hr"
    )
    ax.grid(True)
    ax.set_xlim(0.0, 1.0)
    ax.set_ylim(100.0, 300.0)
    ax.legend(loc="center left", bbox_to_anchor=(1.02, 0.5))
    fig.tight_layout()

    outpng = f"plot/{scheme}_{tag}.png"
    fig.savefig(outpng, dpi=200, bbox_inches="tight")
    plt.close(fig)


def plot_error_compare(dt, tag, t_target):
    """
    Mirrors compare_error_schemes.gp

    Original gnuplot:
      idx = round(t_target / 0.1)
    """
    idx = round(t_target / 0.1)

    datasets = [
        ("FTCS explicit",  f"data/error_ftcs_explicit_{tag}.txt", "o"),
        ("FTCS implicit",  f"data/error_ftcs_implicit_{tag}.txt", "s"),
        ("Dufort Frankel", f"data/error_dufort_{tag}.txt", "^"),
        ("Crank Nicolson", f"data/error_cn_{tag}.txt", "d"),
    ]

    fig, ax = plt.subplots(figsize=(10, 7))

    for label, fname, marker in datasets:
        blocks = read_blocks(fname)
        if idx < len(blocks):
            b = blocks[idx]
            ax.plot(
                b["x"], b["y"],
                marker=marker,
                linewidth=1.2,
                markersize=6,
                label=label
            )

    # horizontal zero line, like gnuplot's "0 with lines"
    ax.axhline(0.0, linewidth=2, linestyle="--", label="Exact error (0)")

    ax.set_xlabel("X [ft]")
    ax.set_ylabel("Error = T_scheme - T_exact [°F]")
    ax.set_title(f"Error vs Exact, dt = {dt:.2f} hr, t = {t_target:.1f} hr")
    ax.grid(True)
    ax.set_xlim(0.0, 1.0)
    ax.legend(loc="center left", bbox_to_anchor=(1.02, 0.5))
    fig.tight_layout()

    outpng = f"plot/error_schemes_{tag}_t{t_target:.1f}.png"
    fig.savefig(outpng, dpi=200, bbox_inches="tight")
    plt.close(fig)


def read_convergence_file(filename):
    data = np.loadtxt(filename, comments="#")
    return {
        "dt": data[:, 0],
        "L2_exp": data[:, 1],
        "L2_imp": data[:, 2],
        "L2_df": data[:, 3],
        "L2_cn": data[:, 4],
        "Linf_exp": data[:, 5],
        "Linf_imp": data[:, 6],
        "Linf_df": data[:, 7],
        "Linf_cn": data[:, 8],
    }


def plot_convergence(t_target):
    """
    Mirrors convergence.gp
    Only plots L2 columns, same as your active gnuplot script.
    """
    ttag = make_t_tag(t_target)
    infile = f"data/convergence_t{ttag}.txt"
    conv = read_convergence_file(infile)

    fig, ax = plt.subplots(figsize=(10, 7))

    ax.loglog(conv["dt"], conv["L2_exp"], marker="o", linewidth=2, markersize=6, label="FTCS explicit L2")
    ax.loglog(conv["dt"], conv["L2_imp"], marker="s", linewidth=2, markersize=6, label="FTCS implicit L2")
    ax.loglog(conv["dt"], conv["L2_df"], marker="^", linewidth=2, markersize=6, label="Dufort Frankel L2")
    ax.loglog(conv["dt"], conv["L2_cn"], marker="d", linewidth=2, markersize=6, label="Crank Nicolson L2")

    ax.set_xlabel("dt [hr]")
    ax.set_ylabel("Error norm")
    ax.set_title(f"Convergence vs dt at t = {t_target:.2f} hr")
    ax.grid(True, which="both")
    ax.legend(loc="center left", bbox_to_anchor=(1.02, 0.5))
    fig.tight_layout()

    outpng = f"plot/convergence_t{t_target:.2f}.png"
    fig.savefig(outpng, dpi=200, bbox_inches="tight")
    plt.close(fig)


def main():
    if len(sys.argv) != 3:
        print("Usage: python3 plot_heat_1d.py <dt> <t_target>")
        sys.exit(1)

    dt = float(sys.argv[1])
    t_target = float(sys.argv[2])

    tag = make_dt_tag(dt)

    ensure_plot_dir()

    plot_scheme(dt, "ftcs_explicit", tag)
    plot_scheme(dt, "dufort", tag)
    plot_scheme(dt, "ftcs_implicit", tag)
    plot_scheme(dt, "cn", tag)

    plot_error_compare(dt, tag, t_target)
    plot_convergence(t_target)

    print("Plots created successfully.")


if __name__ == "__main__":
    main()