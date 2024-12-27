import numpy as np
import matplotlib.colorbar as colorbar
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
import matplotlib.cm as cm
import os


def tanh_func(x, A, B, C, D):
    """Hyperbolic tangent function for curve fitting."""
    return A * np.tanh(B * (x - D)) + C

def analyse_profile(filename):
    data = np.loadtxt(filename)
    symmetrized_data = 0.5 * (data + data[::-1])
    half_data = symmetrized_data[:25]

    popt, _ = curve_fit(tanh_func, np.arange(25),half_data, p0=[0.1, 0.1, 0.1, 12])
    A, B, C, D = popt
    return 2*D - 2*1.83/B


print(analyse_profile('volume_frac_pol_9.4_10Beads.txt'), ' boundary center and width for compact spacers')

print(analyse_profile('volume_frac_pol_9_450Beads.txt'), ' boundary center and width for IDR spacers')
