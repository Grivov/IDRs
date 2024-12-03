import numpy as np
import matplotlib.pyplot as plt
import os


def load_and_symmetrize(filename):
    """ Load data from a file and symmetrize it. """
    data = np.loadtxt(filename)
    # Symmetrize the data
    symmetrized_data = 0.5 * (data + data[::-1])
    return symmetrized_data

def main():
    # Base directory where files are stored
    base_dir = os.path.expanduser("~/Desktop/volume_fractions/")
    

    # File names
    filenames = [
        "volume_frac_pol_9.4_10Beads.txt",
        "volume_frac_pol_9_450Beads.txt",
        "volume_frac_par_9.4_10Beads.txt",
        "volume_frac_par_9_450Beads.txt"
    ]    
    
    # Plot titles
    titles = [
        "volume fraction of polymers with large spacers",
        "volume fraction of polymers with small spacers",
        "volume fraction of intruders",
        "volume fraction of intruders"
    ]
    
    # Y-axis limits for each plot
    y_limits = [
        (0, 0.1),  # Upper left
        (0, 0.1),  # Upper right
        (0, 0.019),  # Lower left
        (0, 0.025)   # Lower right
    ]
    
    # Setting plot aesthetics
    plt.rcParams['font.size'] = 20  # Increase font size for all text
    plt.rcParams['axes.titlesize'] = 20  # Font size for titles
    plt.rcParams['axes.labelsize'] = 20  # Font size for labels
    plt.rcParams['xtick.labelsize'] = 20  # Font size for X tick labels
    plt.rcParams['ytick.labelsize'] = 20  # Font size for Y tick labels
    plt.rcParams['lines.linewidth'] = 5  # Increase line width
    
    # Create a 2x2 subplot grid
    fig, axs = plt.subplots(2, 2, figsize=(40, 10))
    plt.subplots_adjust(wspace=0.15, hspace=0.4)
    x_values = np.linspace(-50, 50, 50)
    # Iterate over files and plot each in its respective subplot
    for ax, filename, title, ylim in zip(axs.ravel(), filenames, titles, y_limits):
        full_path = os.path.join(base_dir, filename)
        data = load_and_symmetrize(full_path)
        
        # Determine plot color and style based on position
        if ylim[1] > 0.025:
            color = 'black'
        else:
            #color = '#CCFF66'
            color = '#74B72E'
        
        # Plotting the data with both line and scatter plot
        ax.plot(x_values,data, color=color, linewidth=7)  # Line plot
        ax.scatter(x_values, data, color=color, s=200)  # Scatter plot
        
        #ax.set_title(title)
        ax.set_ylim(ylim)  # Set y-axis limits
        
    #plt.tight_layout()
    plt.savefig('volume_fractions_green.png', dpi = 800)
    plt.show()

if __name__ == "__main__":
    main()