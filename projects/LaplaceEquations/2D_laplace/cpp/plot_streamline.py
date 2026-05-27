import matplotlib.pyplot as plt
import numpy as np

data = np.loadtxt("data/psi.dat")

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
plt.title("Streamlines from stream function")
plt.axis("equal")
plt.savefig("plot/streamline.png", dpi=300)
