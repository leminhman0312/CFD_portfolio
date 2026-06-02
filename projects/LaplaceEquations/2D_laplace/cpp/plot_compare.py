import matplotlib.pyplot as plt
import numpy as np

solvers = [
    ("point_gs", "Point GS"),
    ("point_sor", "Point SOR"),
    ("line_gs", "Line GS"),
    ("line_sor", "Line SOR"),
]

plt.figure()

for filename, label in solvers:
    data = np.loadtxt(f"data/{filename}_history.dat")

    iteration = data[:, 0]
    error = data[:, 1]

    plt.semilogy(iteration, error, linewidth=2, label=label)

plt.xlabel("Iteration")
plt.ylabel("Max error")
plt.title("Convergence history comparison")

plt.grid(True)
plt.legend()

plt.savefig("plot/convergence_comparison.png", dpi=300, bbox_inches="tight")

print("Saved plot/convergence_comparison.png")
