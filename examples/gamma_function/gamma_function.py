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

def plot_gamma_function():
    # File path
    file_path = os.path.join(os.path.dirname(__file__), 'gamma_function.dat')
    
    if not os.path.exists(file_path):
        print(f"Error: {file_path} not found.")
        return

    # Load data
    # Header: # x Gamma(x)
    data = np.loadtxt(file_path)
    x = data[:, 0]
    y = data[:, 1]

    # Create figure
    fig, ax = plt.subplots(figsize=(10, 8), constrained_layout=True)
    fig.suptitle(r'Gamma Function $\Gamma(x)$ Analysis', fontsize=20, fontweight='bold', color='#00d4ff', y=1.05)

    # Plot the function
    # To avoid connecting lines across poles, we split by jumps in x or use a threshold
    # Since we have segments (between negative integers), let's detect them
    dx = np.diff(x)
    median_dx = np.median(dx)
    split_indices = np.where(dx > median_dx * 5)[0] + 1
    
    x_segments = np.split(x, split_indices)
    y_segments = np.split(y, split_indices)

    for i, (xs, ys) in enumerate(zip(x_segments, y_segments)):
        label = r'$\Gamma(x)$' if i == 0 else None
        ax.plot(xs, ys, color='#ff007c', linewidth=2.5, label=label)

    # Aesthetics
    ax.set_title('Real Gamma Function Visualization', fontsize=14, pad=10)
    ax.set_xlabel('x', fontsize=12)
    ax.set_ylabel(r'$\Gamma(x)$', fontsize=12)
    ax.grid(True, alpha=0.3, linestyle='--')
    
    # Limits and reference lines
    ax.set_ylim(-10, 10)
    ax.set_xlim(min(x), max(x))
    ax.axhline(0, color='white', linewidth=0.8, alpha=0.5)
    
    # Highlight poles
    poles = [0, -1, -2, -3, -4]
    for pole in poles:
        if min(x) <= pole <= max(x):
            ax.axvline(pole, color='#00ff9d', linestyle=':', alpha=0.5, linewidth=1.5)

    ax.legend()

    # Save the figure
    output_path = os.path.join(os.path.dirname(__file__), 'gamma_function_plot.png')
    plt.savefig(output_path, dpi=300, bbox_inches='tight')
    print(f"Plot saved to {output_path}")

if __name__ == "__main__":
    plot_gamma_function()
