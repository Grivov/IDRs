import numpy as np
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker

# [Previous data definitions remain the same]
partitioning_large = [
    9.357753842811536993e-01,
    7.975284462409490960e-01,
    5.967281665983523720e-01,
    3.854557569297885622e-01,
    2.151985116763894257e-01,
    1.056790277684135637e-01,
    4.664814051407292078e-02,]

partitioning_small = [
    9.296126177707715321e-01,
    5.460457766262453161e-01,
    2.258494605550781520e-01,
    7.175230556080981759e-02,
    1.837496318223453065e-02,
    3.943210090582702741e-03,
    7.667199030448663315e-04]

partitioning_box = [
    9.298590493621123887e-01,
    5.123187132264528376e-01,
    1.485373478812549441e-01,
    2.054308042807257914e-02,
    1.215484434864740320e-03,
    2.728102990929338019e-05,
    2.231193511229733852e-07]

partitioning_box_2 = [partitioning*np.exp((4 * np.pi * (R)**3 / 3) * 0.17/ 4.14) 
                      for R, partitioning in zip(np.arange(0, len(partitioning_small) * 0.5, 0.5), partitioning_box)]

pressure_box = 0.1715
green = '#77A632'
orng = "#D9501E"

n = np.arange(0, len(partitioning_small) * 0.5, 0.5)
log_partitioning_large = -np.log(partitioning_large)
log_partitioning_small = -np.log(partitioning_small)
log_partitioning_box_2 = -np.log(partitioning_box_2)

# Create the plot with subplots
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(16, 6), gridspec_kw={'wspace': 0.25})

# Plot for the logarithmic transformation
ax1.plot(n[1:], log_partitioning_large[1:], color=green, marker='s', markersize=22, linewidth=5, label='Compact spacers')
ax1.plot(n[1:], log_partitioning_small[1:], color=orng, marker='o', markersize=22, linewidth=5, label='IDR spacers', zorder = 9)
#ax1.plot(n[1:], log_partitioning_box_2[1:], color=orng, marker='o', markersize=22, linewidth=5, label='Spacers-only')
ax1.scatter(n[1:], log_partitioning_box_2[1:], color='white', marker='o', s=100, linewidth=5, zorder = 8)
# Configure axis 1
ax1.set_xlabel('R, nm', fontsize=19)
ax1.set_ylabel('∆F, kT', fontsize=19)
ax1.set_xscale('log')
ax1.set_yscale('log')

# Set custom ticks
ytick_values = np.exp(np.linspace(np.log(0.15), np.log(10), 7))
xtick_values = np.linspace(0.5, 2.5, 5)

# Format tick labels
ytick_labels = [f'{val:.1f}' for val in ytick_values]
xtick_labels = [f'{val:.1f}' for val in xtick_values]

# Apply ticks
ax1.set_yticks(ytick_values)
ax1.set_yticklabels(ytick_labels)
ax1.set_xticks(xtick_values)
ax1.set_xticklabels(xtick_labels)

# Remove minor ticks
ax1.minorticks_off()

# Set other properties
ax1.tick_params(axis='both', labelsize=19)
ax1.legend(fontsize=19)

# Plot for original data
ax2.plot(n, partitioning_large, color=green, marker='s', markersize=22, linewidth=5, label='Compact spacers')
ax2.plot(n, partitioning_small, color=orng, marker='o', markersize=22, linewidth=5, label='IDR spacers', zorder = 9)
#ax2.plot(n, partitioning_box_2, color=orng, marker='o', markersize=22, linewidth=5, label='Spacers-only')
ax2.scatter(n, partitioning_box_2, color='white', marker='o', s=100, linewidth=5, zorder = 8)

# Add large plus markers
ax2.plot(1.5, 0.0736, color='black', marker='+', markersize=40, markeredgewidth=3, linestyle='None', zorder = 10)
ax2.plot(1.5, 0.386, color='black', marker='+', markersize=40, markeredgewidth=3, linestyle='None', zorder = 10)

# Configure axis 2
ax2.set_xlabel('R, nm', fontsize=19)
ax2.set_ylabel('P(R)', fontsize=19)
ax2.tick_params(axis='both', labelsize=19)
ax2.legend(fontsize=19)

plt.tight_layout()
plt.savefig('small_large_spacers.png', dpi=600, bbox_inches='tight')
plt.show()