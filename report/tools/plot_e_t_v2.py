import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
import matplotlib.patches as patches
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('-i', '--input', type=str, help='input file', required=True)
# parser.add_argument('-o', '--output', type=str, help='output file', required=True)
parser.add_argument('-s', '--highlight_start', type=float, help='highlight start', required=True)
parser.add_argument('-e', '--highlight_end', type=float, help='highlight end', required=True)
args = parser.parse_args()

# Plt settings
plt.rcParams['font.family'] = 'Arial'
sns.set_style("whitegrid")

# Get data
data = np.loadtxt(args.input)
energy = data[:, 0]
temperature = data[:, 1]

# Generate sample data
# energy = np.linspace(0, 10, 100)  # Total energy values in eV
# temperature = 300 + np.random.randn(100) * 50  # Temperature values in K

# Create the main plot
plt.figure(figsize=(10, 6))
sns.lineplot(x=energy, y=temperature, markers=True)
plt.xlabel('Total Energy (eV)')
plt.ylabel('Temperature (K)')
plt.title('Temperature vs. Total Energy')

# Highlight a certain range on the main plot
highlight_start = args.highlight_start
highlight_end = args.highlight_end
temperature_start = temperature[np.where(energy < highlight_start)[0][-1]]
temperature_end = temperature[np.where(energy > highlight_end)[0][0]]
plt.axvspan(highlight_start, highlight_end, color='gray', alpha=0.3)

# # Create a gray zoom-in shape
# rect = patches.Rectangle(
#     (highlight_start, temperature_start),   # (x,y) starting point of the rectangle
#     highlight_end - highlight_start,          # width
#     temperature_end - temperature_start, # height
#     linewidth=1, edgecolor='gray', facecolor='gray', alpha=0.2
# )
# plt.gca().add_patch(rect)
# 
# Create a zoomed-in subplot
zoom_ax = plt.axes([0.8, 0.2, 0.1, 0.1]) # left, bottom, width, height  # Adjust position and size as needed
sns.lineplot(x=energy, y=temperature, ax=zoom_ax)
zoom_ax.set_xlim(highlight_start, highlight_end)
# Limit y axis to the range of temperatures in the zoomed-in area
# Get first data point before highlight_start
# Get temperature
zoom_ax.set_ylim(temperature_start, temperature_end)
zoom_ax.set_xlabel('Total Energy (eV)')
zoom_ax.set_ylabel('Temperature (K)')
zoom_ax.set_title('Phase Transition')

# Display the plot
plt.tight_layout()
plt.show()
