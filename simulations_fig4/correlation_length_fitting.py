import numpy as np
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt
from matplotlib.ticker import MaxNLocator


press = np.loadtxt('/home/vgrigor2/scratch4-yzhan567/vgrigor/idr/figure3_boxes/multiple_boxes_3to8/pressures.txt')

nu = 0.588

press = np.array(press)
xi = (4.14 / (3 * nu - 1) / press) ** (1/3)
print(xi)

vf = np.exp(np.linspace(np.log(0.03), np.log(0.08),10))

# Fit the data to the form A*vf**B
def func(vf, A, B):
    return A * vf**B

popt, pcov = curve_fit(func, vf, xi)
A, B = popt
xi_fit = func(vf, A, B)

# Plot the data points
plt.figure(figsize=(12, 10))
plt.plot(vf, xi, 'o', label='Data', markersize = 20)


def saw(vf, c):
    return c * vf ** (-nu / (3 * nu - 1))

popt2, pcov2 = curve_fit(saw, vf, xi)
c = popt2

# Plot the line
#c = 1.06*0.38 # Adjust this value to change the scale of the line
#xi_line = c * vf ** (-0.5 / (3 * 0.5 - 1))

vf_fine = np.linspace(0.023,0.08,100)
xi_line = saw(vf_fine, c)

plt.plot(vf_fine, xi_line, '-', label='self avoiding walk')

 # Ensure x-axis ticks are integers
ax = plt.gca()  # Get the current Axes instance on the current figure matching the given keyword args, or create one.
ax.xaxis.set_major_locator(MaxNLocator(integer=True))

#B*3*nu-B = -nu 
#nu(3*B+1) = B
#nu = B/(3*B+1)


# Plot the fitted line
plt.loglog(vf, xi_fit, '--', label=f'fitted line')
#plt.plot(vf, xi_fit, '--', label=f'fitted line')

# Add labels and legend
plt.xlabel('phi')
plt.ylabel('xi')
plt.title('Data and Line Plot')
plt.legend()
plt.grid()


plt.show()
