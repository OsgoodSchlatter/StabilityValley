import matplotlib.pyplot as plt
import numpy as np
import matplotlib.colors as mcolors
from matplotlib.patches import Patch

# Define the list of atoms with their data
atoms = [
    [
        [["hydrogen1"], ["H"], [1], [1], [0], [-1], ["st"]],
        [["hydrogen2"], ["H"], [1], [2], [1112], [-1], ["st"]],
        [["hydrogen3"], ["H"], [1], [3], [2827], [388781328], ["b-"]],
        [["hydrogen4"], ["H"], [1], [4], [1720], [0], ["n"]],
        [["hydrogen5"], ["H"], [1], [5], [1336], [86 * 10**(-24)], ["n"]],
        [["hydrogen6"], ["H"], [1], [6], [960], [290 * 10**(-24)], ["none"]],
        [["hydrogen7"], ["H"], [1], [7], [940], [10 * 10**(-21)], ["none"]],
    ],
    [
        [["helium3"], ["He"], [2], [3], [2572], [-1], ["st"]],
        [["helium4"], ["He"], [2], [4], [7074], [-1], ["st"]],
        [["helium5"], ["He"], [2], [5], [5512], [704*10 ** (-24)], ["n"]],
        [["helium6"], ["He"], [2], [6], [4878], [0.8], ["b-"]],
        [["helium7"], ["He"], [2], [7], [4123], [3 * 10**(-21)], ["n"]],
        [["helium8"], ["He"], [2], [8], [3924], [0.119], ["b-"]],
        [["helium9"], ["He"], [2], [9], [3349], [0], ["n"]],
        [["helium10"], ["He"], [2], [10], [2995], [1.5*10**(-21)], ["n"]]
    ]
]

# Create the figure and axes
fig1, ax1 = plt.subplots(figsize=(10, 6))
fig2, ax2 = plt.subplots(figsize=(10, 6))
fig3, ax3 = plt.subplots(figsize=(10, 6))

# Define a custom colormap for the scale
scale_colormap = plt.cm.viridis

# Initialize variables to track the minimum and maximum values for half-life and binding energy
min_half_life = float('inf')
max_half_life = float('-inf')
min_binding_energy = float('inf')
max_binding_energy = float('-inf')

# Iterate over each set of atoms
for atom_data in atoms:
    # Extract data for the current set of atoms
    A_minus_Z = [item[3][0] - item[2][0] for item in atom_data]
    Z_values = [item[2][0] for item in atom_data]
    binding_energy = [item[4][0] for item in atom_data]
    half_life = [item[5][0] for item in atom_data]
    decay_types = [item[6][0] for item in atom_data]

    # Update min and max values for half-life and binding energy
    min_half_life = min(min_half_life, min(half_life))
    max_half_life = max(max_half_life, max(half_life))
    min_binding_energy = min(min_binding_energy, min(binding_energy))
    max_binding_energy = max(max_binding_energy, max(binding_energy))

    # Define colors based on binding energy, half-life, or decay type
    color_map_binding_energy = plt.cm.viridis(
        (np.array(binding_energy) - min_binding_energy) / (max_binding_energy - min_binding_energy))

    half_life_colors = []
    for value in half_life:
        if value < 0:
            half_life_colors.append((0, 0, 0, 1))  # Black color
        else:
            half_life_colors.append(plt.cm.plasma(
                (value - min_half_life) / (max_half_life - min_half_life)))

    half_life_colormap = mcolors.ListedColormap(half_life_colors)

    color_map_decay_type = {
        "st": "black",
        "b-": "cyan",
        "n": "pink",
        "none": "gray",
    }

    # Plot data points on the same graph

    # HALF LIFE
    normalized_color_map_hl = plt.Normalize(vmin=0, vmax=max_half_life)
    scatter_halflife = ax1.scatter(A_minus_Z, Z_values, c=half_life,
                                   cmap=scale_colormap, norm=normalized_color_map_hl, marker='s', s=1000, label=' '.join(atom_data[0][1]))

    ax1.set_xlabel("N")
    ax1.set_ylabel("Z")
    ax1.set_title("Half-life of Isotopes")

    # BINDING ENERGY
    normalized_color_map_binding_energy = plt.Normalize(
        vmin=0, vmax=max_binding_energy)

    scatter_energy = ax2.scatter(A_minus_Z, Z_values, c=binding_energy,
                                 cmap=scale_colormap, norm=normalized_color_map_binding_energy, marker='s', s=1000, label=' '.join(atom_data[0][1]))

    ax2.set_xlabel("N")
    ax2.set_ylabel("Z")
    ax2.set_title("Binding Energy of Isotopes")

    # DECAY TYPE
    # For the decay graph, use color_map_decay_type for colors and create a custom legend
    scatter_decay = ax3.scatter(A_minus_Z, Z_values, c=[color_map_decay_type[d] for d in decay_types],
                                cmap=scale_colormap, marker='s', s=1000, label=' '.join(atom_data[0][1]), vmin=0, vmax=1)

    ax3.set_xlabel("N")
    ax3.set_ylabel("Z")
    ax3.set_title("Decay Type of Isotopes")

cbar1 = plt.colorbar(scatter_halflife, ax=ax1,
                     extend='max', label='Half life (s)')
cbar2 = plt.colorbar(scatter_energy, ax=ax2,
                     extend='max', label='Binding Energy (keV)')

# Create a custom legend for the decay graph
legend_labels = list(color_map_decay_type.keys())
legend_elements = [Patch(facecolor=color_map_decay_type[label],
                         edgecolor='black', label=label) for label in legend_labels]
ax3.legend(handles=legend_elements, loc='upper right', title='Decay Types')

# Set the Y-axis limits for all graphs based on Z values
max_Z = max(max(Z_values) for atom_data in atoms)
ax1.set_ylim(0, max_Z + 1)
ax2.set_ylim(0, max_Z + 1)
ax3.set_ylim(0, max_Z + 1)

# Show the plots
plt.show()
