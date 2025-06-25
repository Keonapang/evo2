
import matplotlib.pyplot as plt
import numpy as np
from mpl_toolkits.mplot3d import Axes3D

# Creating the values for x and y axis
x = [1, 2, 3, 4, 5, 6, 7]
y = [1, 2, 3, 4]
x, y = np.meshgrid(x, y)

# AUC values
t1024 = np.array([0.8425, 0.886948529, 0.879032258, 0.7328, 0.744653518, 0.743377282, 0.739842661])
t2048 = np.array([0.945, 0.953125, 0.913151365, 0.8235, 0.824025794, 0.806951729, 0.804288383])
t4096 = np.array([0.9525, 0.953125, 0.910049628, 0.8346, 0.836570191, 0.811037053, 0.806454966])
t8192 = np.array([0.9525, 0.953125, 0.925248139, 0.8739, 0.873321746, 0.870751719, 0.866670692])
z = np.array([t1024, t2048, t4096, t8192])

# Create the 3D plot
fig = plt.figure(figsize=(25, 30))  # Increase figure size for more white space
ax = plt.axes(projection='3d')

# Use a surface plot for better visualization
surf = ax.plot_surface(x, y, z, cmap='viridis', edgecolor='none')

# Customize ticks and labels
new_tick_locations = [1, 2, 3, 4, 5, 6, 7]
new_tick_labels = ['50', '84', '150', '250', '500', '1000', '2000']
ax.set_xticks(new_tick_locations)
ax.set_xticklabels(new_tick_labels, fontsize=15)  # Increase font size of tick labels

new_tick_locations2 = [1, 2, 3, 4]
new_tick_labels2 = ['1024', '2048', '4096', '8192']
ax.set_yticks(new_tick_locations2)
ax.set_yticklabels(new_tick_labels2, fontsize=15)

# Set labels and title
ax.set_xlabel('Input seq length', labelpad=25, fontsize=20)  # Add padding to move the label farther
ax.set_ylabel('Window Size', labelpad=25, fontsize=20)      # Add padding to move the label farther
ax.set_zlabel('AUC', labelpad=25, fontsize=20)              # Add padding to move the label farther
ax.set_title('AUC Performance of Evo2 7B Model', fontweight='bold', pad=30, fontsize=20)  # Bold title and padding

ax.tick_params(axis='z', labelsize=14)  # Make z-axis tick labels larger
ax.tick_params(axis='x', pad=11)        # Move x-axis tick labels farther from the plot
ax.tick_params(axis='y', pad=8)        # Move x-axis tick labels farther from the plot
ax.tick_params(axis='z', pad=11)        # Move x-axis tick labels farther from the plot

# Adjust viewing angle
ax.view_init(50, 70)
# Add a color bar to indicate AUC values with reduced size
cbar = fig.colorbar(surf, ax=ax, shrink=0.4, aspect=10)  # Shrink and reduce aspect ratio
cbar.ax.tick_params(labelsize=15)  # Smaller font size for color bar ticks

# Adjust margins to add more white space
plt.subplots_adjust(left=0.25, right=0.8, top=0.85, bottom=0.2)

# Save the plot
output_file = "/mnt/nfs/rigenenfs/workspace/pangk/Softwares/evo2/results/auc_3d_plot_adjusted.png"
plt.savefig(output_file, dpi=300, bbox_inches='tight')  # Save with high resolution
plt.show()
print(f"The plot has been saved to {output_file}")

