import sys

import matplotlib.pyplot as plt
import numpy as np

solver = sys.argv[1]

datafile = f"data/{solver}.dat"
plotfile = f"plot/{solver}.png"

data = np.loadtxt(datafile)

x = data[:, 0]
y = data[:, 1]
psi = data[:, 2]

imax = 31
jmax = 21

X = x.reshape(jmax, imax)
Y = y.reshape(jmax, imax)
PSI = psi.reshape(jmax, imax)

plt.figure()
plt.contour(X, Y, PSI, levels=20)
plt.xlabel("x")
plt.ylabel("y")
plt.title(solver)
plt.axis("equal")

plt.savefig(plotfile, dpi=300)
print(f"Saved {plotfile}")
