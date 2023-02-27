import os
import argparse
import logging
from features import DataSetMLP

logging.getLogger("tensorflow").disabled=True
logging.getLogger("h5py._conv").disabled=True

# read data and convert to binary format

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-f", "--file", type=str, help="path to the data file")
    args = parser.parse_args()

    work_directory = f"{os.path.dirname(args.file)}"

    data = DataSetMLP(file=args.file)
    data.save_inputs(f"{work_directory}/inputs.npy")
    data.save_labels(f"{work_directory}/labels.npy")