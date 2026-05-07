#!/usr/bin/env python3

import os
import sys

import matplotlib.pyplot as plt
import numpy as np


def read_xyz(filename):
    data = np.loadtxt(filename)

    x = data[:, 0]
    y = data[:, 1]
    T = data[:, 2]

    nx = len(np.unique(x))
    ny = len(np.unique(y))

    X = x.reshape(ny, nx)
    Y = y.reshape(ny, nx)
    Z = T.reshape(ny, nx)

    return X, Y, Z


def plot_contour(datafile, outpng, time_hr, scheme):
    X, Y, Z = read_xyz(datafile)

    os.makedirs(os.path.dirname(outpng), exist_ok=True)

    fig, ax = plt.subplots(figsize=(8, 7))

    cf = ax.contourf(X, Y, Z, levels=40)
    ax.contour(X, Y, Z, levels=12, linewidths=0.5)

    cbar = fig.colorbar(cf, ax=ax)
    cbar.set_label("Temperature")

    ax.set_xlabel("x")
    ax.set_ylabel("y")
    ax.set_title(f"{scheme}, t = {time_hr:.3f}")
    ax.set_aspect("equal", adjustable="box")

    fig.tight_layout()
    fig.savefig(outpng, dpi=200)
    plt.close(fig)


def plot_convergence(datafile, outpng):
    data = np.loadtxt(datafile)

    dt = data[:, 0]
    l2 = data[:, 1]
    linf = data[:, 2]

    os.makedirs(os.path.dirname(outpng), exist_ok=True)

    fig, ax = plt.subplots(figsize=(8, 6))

    ax.loglog(dt, l2, marker="o", linewidth=2, label="L2 error")
    ax.loglog(dt, linf, marker="s", linewidth=2, label="Linf error")

    ax.set_xlabel("dt")
    ax.set_ylabel("error")
    ax.set_title("Time convergence, implicit ADI")
    ax.grid(True, which="both")
    ax.legend()

    fig.tight_layout()
    fig.savefig(outpng, dpi=200)
    plt.close(fig)


def main():
    if len(sys.argv) < 4:
        print("Usage:")
        print("  python3 plot_heat_2d.py contour <datafile> <outpng> <time> <scheme>")
        print("  python3 plot_heat_2d.py convergence <datafile> <outpng>")
        sys.exit(1)

    mode = sys.argv[1]

    if mode == "contour":
        plot_contour(sys.argv[2], sys.argv[3], float(sys.argv[4]), sys.argv[5])
    elif mode == "convergence":
        plot_convergence(sys.argv[2], sys.argv[3])
    else:
        print(f"Unknown mode: {mode}")
        sys.exit(1)


if __name__ == "__main__":
    main()
