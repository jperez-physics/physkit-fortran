import os
import numpy as np
import matplotlib.pyplot as plt

# ── Data ──────────────────────────────────────
base_path = os.path.dirname(os.path.abspath(__file__))
file = os.path.join(base_path, "gamma_graph.dat")

col_x   = 0
col_y   = 1
# ─────────────────────────────────────────────

data = np.loadtxt(file, comments="#")

plt.plot(data[:, col_x], data[:, col_y])
plt.ylim(-2, 7)
plt.xlabel("X")
plt.ylabel("Y")
plt.title(os.path.basename(file))
plt.grid(True)
plt.tight_layout()
plt.savefig('gamma_graph.png')
plt.show()