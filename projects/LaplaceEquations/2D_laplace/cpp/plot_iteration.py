import numpy as np
import matplotlib.pyplot as plt

data = np.loadtxt("data/line_sor_w_sweep.dat")

w = data[:,0]
iters = data[:,1]

plt.figure()
plt.plot(w, iters, marker="o")

plt.xlabel("Relaxation factor, w")
plt.ylabel("Iterations")
plt.title("Line SOR: Iterations vs Relaxation Factor")

plt.grid(True)

plt.savefig("plot/line_sor_w_sweep.png", dpi=300)
