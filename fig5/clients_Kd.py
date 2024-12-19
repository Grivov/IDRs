import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import minimize
from scipy.optimize import curve_fit
from scipy.stats import chi2
import matplotlib.cm as cm
from scipy import integrate
from scipy.interpolate import CubicSpline
plt.rcParams.update({'font.size': 20})

P_values = np.array([0.077, 0.080, 0.077, 0.083, 0.102, 0.147, 0.521, 2.011, 16.107, 2384.184])

Kd = np.array([3.252331479106505441e+00,
1.964557011599213254e+00,
1.165642225778261709e+00,
6.449599101842320437e-01,
2.517266024767562671e-01,
7.519506129417723173e-02,
2.665710453967117671e-02,
1.590682305989360099e-02,
1.116392790971670879e-02,
7.371598581233341878e-03,])

rho_array = np.array([6.858732348524217148e-02,
6.845226617284526360e-02,
6.879947066181733850e-02,
6.831391479004834955e-02,
6.740303499215570537e-02,
6.786648036744238888e-02,
6.504726408128631843e-02,
6.369209432737790721e-02,
6.622419709877211402e-02,
7.424716202761622030e-02,])

F_req = - rho_array / Kd

plt.figure(figsize=(20, 5))
cmap = cm.get_cmap('viridis', 10)
colors_indices = np.arange(10)

# Create scatter plots with different labels for each epsilon
for i in range(len(colors_indices)):
    plt.scatter(Kd[i]*0.6, -np.log(P_values[i]), c=[cmap(i/9)], s=600, zorder=5, 
               label=f'ε = {i+1}  ')

plt.plot(Kd*0.6, F_req -np.log(7.17e-02),  linewidth=5, c='black', alpha=0.7)
plt.scatter(Kd*0.6, F_req -np.log(7.17e-02), label='prediction', linewidth=5, c='black', s=400)
plt.xscale('log')

# Add legend in lower right
plt.legend(loc='lower right', fontsize=25, ncol=3)
plt.savefig('clients_Kd_curve', dpi = 500)

plt.show()
