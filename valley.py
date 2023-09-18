import matplotlib.pyplot as plt
import numpy as np
import matplotlib.colors as mcolors
from matplotlib.patches import Patch

# Your hydrogen data
hydrogen = [
    [["hydrogen1"], ["H"], [1], [1], [0], [-1], ["st"]],
    [["hydrogen2"], ["H"], [1], [2], [1112], [-1], ["st"]],
    [["hydrogen3"], ["H"], [1], [3], [2827], [388781328], ["b-"]],
    [["hydrogen4"], ["H"], [1], [4], [1720], [0], ["n"]],
    [["hydrogen5"], ["H"], [1], [5], [1336], [86 * 10**(-24)], ["n"]],
    [["hydrogen6"], ["H"], [1], [6], [960], [290 * 10**(-24)], ["none"]],
    [["hydrogen7"], ["H"], [1], [7], [940], [10 * 10**(-21)], ["none"]],
]

# Extract data
A_minus_Z = [item[3][0] - item[2][0] for item in hydrogen]
Z_values = [item[2][0] for item in hydrogen]
binding_energy = [item[4][0] for item in hydrogen]
half_life = [item[5][0] for item in hydrogen]
decay_types = [item[6][0] for item in hydrogen]

# Define colors based on binding energy, half-life, or decay type
color_map_binding_energy = plt.cm.viridis(
    np.array(binding_energy) / max(binding_energy))

half_life_colors = []
for value in half_life:
    if value < 0:
        half_life_colors.append((0, 0, 0, 1))  # Black color
    else:
        half_life_colors.append(plt.cm.plasma(value / max(half_life)))

half_life_colormap = mcolors.ListedColormap(half_life_colors)

color_map_decay_type = {
    "st": "black",
    "b-": "cyan",
    "n": "pink",
    "none": "gray",
}

# Create the figure and axes
fig1, ax1 = plt.subplots(figsize=(10, 6))
fig2, ax2 = plt.subplots(figsize=(10, 6))
fig3, ax3 = plt.subplots(figsize=(10, 6))

# Plot data points and create legends

# HALF LIFE
scatter_halflife = ax1.scatter(A_minus_Z, Z_values,
                               c=half_life_colors, cmap="viridis", marker='s', s=2500)
cbar1 = plt.colorbar(scatter_halflife)
cbar1.set_label("Half life (s)")

ax1.set_xlabel("N")
ax1.set_ylabel("Z")
ax1.set_title("Hydrogen Isotopes")

# BINDING ENERGY
scatter_binding = ax2.scatter(A_minus_Z, Z_values,
                              c=color_map_binding_energy, cmap="viridis", marker='s', s=2500)
cbar2 = plt.colorbar(scatter_binding)
cbar2.set_label("Binding Energy (keV)")

ax2.set_xlabel("N")
ax2.set_ylabel("Z")
ax2.set_title("Hydrogen Isotopes")

# DECAY TYPE
scatter_decay = ax3.scatter(A_minus_Z, Z_values,
                            c=[color_map_decay_type[d] for d in decay_types], cmap="viridis", marker='s', s=2500)

ax3.set_xlabel("N")
ax3.set_ylabel("Z")
ax3.set_title("Hydrogen Isotopes")

# Create custom legend for decay types with square markers filled with colors
legend_labels = list(color_map_decay_type.keys())
legend_colors = [color_map_decay_type[label] for label in legend_labels]
legend_elements = [Patch(facecolor=color_map_decay_type[label],
                         edgecolor='black', label=label) for label in legend_labels]
ax3.legend(handles=legend_elements, loc='upper right', title='Decay Types')

# Show the plots
plt.show()
