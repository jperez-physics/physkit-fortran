import os
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation

# ── Data ──────────────────────────────────────
base_path = os.path.dirname(os.path.abspath(__file__))
file_path = os.path.join(base_path, "three-body-problem.dat")

# Columns: t, x1, y1, x2, y2, x3, y3
data = np.loadtxt(file_path, comments="#")

t = data[:, 0]
x1, y1 = data[:, 1], data[:, 2]
x2, y2 = data[:, 3], data[:, 4]
x3, y3 = data[:, 5], data[:, 6]
# ─────────────────────────────────────────────

# ── Setup Plot ────────────────────────────────
fig, ax = plt.subplots(figsize=(8, 8))
ax.set_aspect('equal')
ax.grid(True, linestyle='--', alpha=0.6)

# Find limits
all_x = np.concatenate([x1, x2, x3])
all_y = np.concatenate([y1, y2, y3])
padding = 0.2
ax.set_xlim(np.min(all_x) - padding, np.max(all_x) + padding)
ax.set_ylim(np.min(all_y) - padding, np.max(all_y) + padding)

ax.set_xlabel("x")
ax.set_ylabel("y")
ax.set_title("Three-Body Problem: Figure-8 Orbit")

# Initialize elements
line1, = ax.plot([], [], 'r-', lw=1, alpha=0.5, label='Body 1')
line2, = ax.plot([], [], 'g-', lw=1, alpha=0.5, label='Body 2')
line3, = ax.plot([], [], 'b-', lw=1, alpha=0.5, label='Body 3')

point1, = ax.plot([], [], 'ro', markersize=8)
point2, = ax.plot([], [], 'go', markersize=8)
point3, = ax.plot([], [], 'bo', markersize=8)

ax.legend(loc='upper right')

# ── Animation ─────────────────────────────────
def init():
    line1.set_data([], [])
    line2.set_data([], [])
    line3.set_data([], [])
    point1.set_data([], [])
    point2.set_data([], [])
    point3.set_data([], [])
    return line1, line2, line3, point1, point2, point3

def update(frame):
    # Historical path (trail)
    line1.set_data(x1[:frame], y1[:frame])
    line2.set_data(x2[:frame], y2[:frame])
    line3.set_data(x3[:frame], y3[:frame])
    
    # Current positions
    point1.set_data([x1[frame]], [y1[frame]])
    point2.set_data([x2[frame]], [y2[frame]])
    point3.set_data([x3[frame]], [y3[frame]])
    
    return line1, line2, line3, point1, point2, point3

# Speed up animation by skipping frames if necessary
step = 2
frames = range(0, len(t), step)

ani = FuncAnimation(fig, update, frames=frames, init_func=init, blit=True, interval=20)

plt.tight_layout()
plt.show()
# ─────────────────────────────────────────────
