import argparse
import os
from data.data import DataSet

'''
    This script reads data from file and converts to binary format
'''

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-f", "--file", type=str, help="path to the data file")
    args = parser.parse_args()
    
    # directory where generated files will be stored
    work_directory = os.path.dirname(args.file)

    # load and save
    data = DataSet(args.file)
    data.save_inputs(f"{work_directory}/inputs")
    data.save_labels(f"{work_directory}/labels")
