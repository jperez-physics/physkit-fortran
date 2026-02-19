import os
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

# ── Data ──────────────────────────────────────
base_path = os.path.dirname(os.path.abspath(__file__))
file = os.path.join(base_path, "gamma_graph.dat")

col_x = 0
col_y = 1
col_z = 2
# ─────────────────────────────────────────────

data = np.loadtxt(file, comments="#")

x = data[:, col_x]
y = data[:, col_y]
z = data[:, col_z]

nx = len(np.unique(x))
ny = len(np.unique(y))

X = x.reshape((nx, ny))
Y = y.reshape((nx, ny))
Z = z.reshape((nx, ny))

# ── Plot 3D ───────────────────────────────────
fig = plt.figure()
ax = fig.add_subplot(projection='3d')

ax.plot_surface(X, Y, Z)

ax.set_xlabel("Re(z)")
ax.set_ylabel("Im(z)")
ax.set_zlabel("Re(Gamma(z))")

plt.title(os.path.basename(file))
plt.tight_layout()
plt.savefig('gamma_surface.png')
plt.show()
