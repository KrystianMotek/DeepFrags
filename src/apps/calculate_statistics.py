import os
import argparse
import numpy as np

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--inputs", type=str)
    parser.add_argument("-l", "--labels", type=str)
    args = parser.parse_args()

    inputs = np.load(args.inputs)
    labels = np.load(args.labels)