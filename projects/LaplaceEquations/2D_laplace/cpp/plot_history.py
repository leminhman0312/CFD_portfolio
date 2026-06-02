import sys

import matplotlib.pyplot as plt
import numpy as np

solver = sys.argv[1]
datafile = f"data/{solver}_history.dat"
plotfile = f"plot/{solver}_history.png"

data = np.loadtxt(datafile)

iteration = data[:, 0]
error = data[:, 1]

plt.figure()

plt.semilogy(iteration, error, linewidth=2)

plt.xlabel("Iteration")
plt.ylabel("Max error")
plt.title(f"{solver} convergence history")

plt.grid(True)

plt.savefig(plotfile, dpi=300, bbox_inches="tight")

print(plotfile)
