import argparse
import os
from features import DataSetMLP

# This script reads data from file and converts to binary format

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-f", "--file", type=str, help="path to the data file")
    args = parser.parse_args()

    name = os.path.splitext(os.path.basename(args.file))[0]

    # directory where generated files will be stored
    work_directory = f"{os.path.dirname(args.file)}/{name}"
    os.makedirs(work_directory)

    data = DataSetMLP(args.file)
    data.save_inputs(f"{work_directory}/inputs.npy")
    data.save_labels(f"{work_directory}/labels.npy")