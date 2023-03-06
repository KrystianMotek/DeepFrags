import argparse
import numpy as np

BOND_LENGTH = 3.8

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-l", "--labels", type=str)
    args = parser.parse_args()
    
    labels = np.load(args.labels)

    # planar and dihedral angles correlation
    # histograms for each secondary structure