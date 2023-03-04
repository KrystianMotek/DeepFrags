import os
import argparse
import numpy as np
from core.features import LabelMLP
from utils.structural import Output

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--inputs", type=str)
    parser.add_argument("-l", "--labels", type=str)
    args = parser.parse_args()

    inputs = np.load(args.inputs)
    labels = np.load(args.labels)

    ss = [LabelMLP.extract_ss(vector) for vector in labels]
    alpha = [Output(vector=vector).alpha() for vector in inputs]
    theta = [Output(vector=vector).theta() for vector in inputs]