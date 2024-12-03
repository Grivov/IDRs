import numpy as np
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt
from matplotlib.ticker import MaxNLocator

press = [2.383746666666666791e-02,
3.055971333333333192e-02,
3.908416666666666983e-02,
5.041458000000000050e-02,
6.520967333333334282e-02,
8.481580000000001063e-02,
1.102572666666666590e-01,
1.433783999999999892e-01,
1.893388000000000015e-01,
2.491207333333333440e-01
]
nu = 0.588
press = np.array(press)
xi = (4.14 / (3 * nu - 1) / press) ** (1/3)
print(xi)

vf = np.exp(np.linspace(np.log(0.03), np.log(0.08), 10))

# Fit the data to the form A*vf**B
def func(vf, A, B):
    return A * vf**B

popt, pcov = curve_fit(func, vf, xi)
A, B = popt
xi_fit = func(vf, A, B)

colors = plt.cm.viridis(np.linspace(0, 1, len(vf)))

# Plot the data points
plt.figure(figsize=(5, 5))
for i, (vf_value, xi_value) in enumerate(zip(vf, xi)):
    plt.plot(vf_value, xi_value, 'o', label='Data', markersize=10, color=colors[i])

def saw(vf, c):
    return c * vf ** (-nu / (3 * nu - 1))

popt2, pcov2 = curve_fit(saw, vf, xi)
c = popt2
print(c, 'C')

vf_fine = np.linspace(0.018, 0.15, 100)
xi_line = saw(vf_fine, c)
plt.loglog(vf_fine, xi_line, '-', c='black')

# Ensure x-axis ticks are integers
ax = plt.gca()
ax.xaxis.set_major_locator(MaxNLocator(integer=True))

plt.xlabel('phi', fontsize=12)
plt.ylabel('xi', fontsize=12)
plt.title('Data and Line Plot')
#plt.legend()
# plt.grid()  # Removed this line
plt.savefig('correlation_length_fitting', dpi=500)
plt.show()