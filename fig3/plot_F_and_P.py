import numpy as np
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
plt.rcParams.update({'font.size': 20})

volume_frac_par_10Beads = np.loadtxt('volume_frac_par_9.4_10Beads.txt')
volume_frac_par_450Beads = np.loadtxt('volume_frac_par_9_450Beads.txt')

volume_frac_par_10Beads = 0.5*(volume_frac_par_10Beads + volume_frac_par_10Beads[::-1])
volume_frac_par_450Beads = 0.5*(volume_frac_par_450Beads + volume_frac_par_450Beads[::-1])

P_sim_large = np.mean(volume_frac_par_10Beads[20:26]) / np.mean(volume_frac_par_10Beads[:6])  
P_sim_small = np.mean(volume_frac_par_450Beads[20:26]) / np.mean(volume_frac_par_450Beads[:6])  

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

green = '#77A632'
orng = "#D9501E"

n = np.arange(0, len(partitioning_small) * 0.5, 0.5)
log_partitioning_large = -np.log(partitioning_large)
log_partitioning_small = -np.log(partitioning_small)

fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(16, 6), gridspec_kw={'wspace': 0.25})

ax1.plot(n, log_partitioning_large, color=green, marker='s', markersize=22, linewidth=5, label='Compact spacers')
ax1.plot(n, log_partitioning_small, color=orng, marker='o', markersize=22, linewidth=5, label='IDR spacers', zorder = 9)
ax1.set_ylim(0, 8)
ax1.legend(fontsize=20)

ax2.plot(n, partitioning_large, color=green, marker='s', markersize=22, linewidth=5, label='Compact spacers')
ax2.plot(n, partitioning_small, color=orng, marker='o', markersize=22, linewidth=5, label='IDR spacers', zorder = 9)

# Add large plus markers
#ax2.plot(1.5, 0.0736, color='black', marker='+', markersize=40, markeredgewidth=3, linestyle='None', zorder = 10)
#ax2.plot(1.5, 0.386, color='black', marker='+', markersize=40, markeredgewidth=3, linestyle='None', zorder = 10)
ax2.plot(1.5, P_sim_large, color='black', marker='+', markersize=40, markeredgewidth=3, linestyle='None', zorder = 10)
ax2.plot(1.5, P_sim_small, color='black', marker='+', markersize=40, markeredgewidth=3, linestyle='None', zorder = 10)


ax2.legend(fontsize=20)

plt.tight_layout()
plt.savefig('small_large_spacers.png', dpi=600, bbox_inches='tight')
plt.show()
