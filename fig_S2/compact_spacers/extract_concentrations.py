import argparse
import numpy as np
import matplotlib.pyplot as plt

def process_volume_fraction(file_name, x):
    """
    Processes a volume fraction data file, symmetrizes it, and calculates averages.
    
    Parameters:
    - file_name: str, path to the input text file containing data.
    - x: int, threshold index for calculating c_dense and c_dilute.
    """
    volume_fraction = np.loadtxt(file_name)
    
    volume_fraction_sym = (volume_fraction + volume_fraction[::-1]) / 2
    
    c_dense = np.mean(volume_fraction_sym[int(25-x):25])   # From index x to 25
    c_dilute = np.mean(volume_fraction_sym[:x])    # From start to index x
    
    # Print results
    print(f"Symmetrized Volume Fraction Data:")
    print(volume_fraction_sym)
    print(f"\nc_dense (mean from index {x} to 25): {c_dense}")
    print(f"c_dilute (mean from start to index {x}): {c_dilute}")
    
    # Plot the original and symmetrized data
    plt.figure(figsize=(8, 6))
    plt.plot(volume_fraction, label="Original Data", linestyle="--", marker="o")
    plt.plot(volume_fraction_sym, label="Symmetrized Data", marker="s")
    plt.axvline(x, color='r', linestyle=':', label=f"x = {x}")
    plt.xlabel("Index")
    plt.ylabel("Volume Fraction")
    plt.title("Symmetrized Volume Fraction Data")
    plt.legend()
    plt.grid()
    #plt.show()

def main():
    parser = argparse.ArgumentParser(description="Process and symmetrize volume fraction data.")
    parser.add_argument("file_name", type=str, help="Path to the input text file containing volume fraction data.")
    parser.add_argument("x", type=int, help="Threshold index for calculating c_dense and c_dilute.")
    args = parser.parse_args()
    process_volume_fraction(args.file_name, args.x)

if __name__ == "__main__":
    main()
