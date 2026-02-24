import numpy as np
import matplotlib.pyplot as plt
import os

plt.style.use('dark_background')
plt.rcParams.update({
    'font.family': 'sans-serif',
    'font.sans-serif': ['Inter', 'Roboto', 'Arial'],
    'axes.facecolor': '#121212',
    'figure.facecolor': '#0a0a0a',
    'axes.edgecolor': '#333333',
    'axes.labelcolor': '#e0e0e0',
    'xtick.color': '#888888',
    'ytick.color': '#888888',
    'grid.color': '#222222',
    'text.color': '#ffffff'
})

def plot_damped_oscillator():
    # File path
    file_path = os.path.join(os.path.dirname(__file__), 'damped_oscillator.dat')
    
    if not os.path.exists(file_path):
        print(f"Error: {file_path} not found.")
        return

    # Load data
    # Header: # Time Position Velocity Acceleration
    data = np.loadtxt(file_path)
    t = data[:, 0]
    x = data[:, 1]
    v = data[:, 2]
    # acc = data[:, 3]

    # Create figure with 3 subplots
    fig, (ax1, ax2, ax3) = plt.subplots(3, 1, figsize=(10, 15), constrained_layout=True)
    fig.suptitle('Damped Harmonic Oscillator Analysis', fontsize=20, fontweight='bold', color='#00d4ff', y=1.02)

    # Plot 1: Position vs Time
    ax1.plot(t, x, color='#ff007c', linewidth=2, label='Position (x)')
    ax1.set_title('Position vs Time', fontsize=14, pad=10)
    ax1.set_xlabel('Time (s)', fontsize=12)
    ax1.set_ylabel('Position (m)', fontsize=12)
    ax1.grid(True, alpha=0.3)
    ax1.legend()

    # Plot 2: Velocity vs Time
    ax2.plot(t, v, color='#00ff9d', linewidth=2, label='Velocity (v)')
    ax2.set_title('Velocity vs Time', fontsize=14, pad=10)
    ax2.set_xlabel('Time (s)', fontsize=12)
    ax2.set_ylabel('Velocity (m/s)', fontsize=12)
    ax2.grid(True, alpha=0.3)
    ax2.legend()

    # Plot 3: Velocity vs Position (Phase Space)
    ax3.plot(x, v, color='#00d4ff', linewidth=2, label='Phase trajectory')
    ax3.set_title('Phase Space (v vs x)', fontsize=14, pad=10)
    ax3.set_xlabel('Position (x)', fontsize=12)
    ax3.set_ylabel('Velocity (v)', fontsize=12)
    ax3.grid(True, alpha=0.3)
    ax3.legend()

    # Save the figure
    plt.savefig(os.path.join(os.path.dirname(__file__), 'damped_oscillator_plots.png'), dpi=300, bbox_inches='tight')
    print(f"Plots saved to {os.path.join(os.path.dirname(__file__), 'damped_oscillator_plots.png')}")

if __name__ == "__main__":
    plot_damped_oscillator()