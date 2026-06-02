import matplotlib.pyplot as plt
import numpy as np

files = [
    "data/centerline_dx04.dat",
    "data/centerline_dx02.dat",
    "data/centerline_dx01.dat",
]

plt.figure()

for filename in files:
    dx = None
    imax = None
    jmax = None

    with open(filename, "r") as f:
        for line in f:
            if line.startswith("# dx"):
                dx = float(line.split()[2])
            elif line.startswith("# imax"):
                imax = int(line.split()[2])
            elif line.startswith("# jmax"):
                jmax = int(line.split()[2])

    data = np.loadtxt(filename)

    x = data[:, 0]
    psi = data[:, 1]

    label = f"dx={dx:g} ({imax}×{jmax})"

    plt.plot(x, psi, marker="o", label=label)

plt.xlabel("x")
plt.ylabel("ψ")
plt.title("Grid refinement study (Line SOR, y = 2.0 m)")

plt.grid(True)
plt.legend()

plt.savefig("plot/grid_refinement.png", dpi=300, bbox_inches="tight")

print("Saved plot/grid_refinement.png")
