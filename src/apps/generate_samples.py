import os
import argparse
from data import all_samples

# get all observations from given directory

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-p", "--path", type=str, help="path to the directory")
    args = parser.parse_args()

    path = args.path

    file = f"{os.path.dirname(path)}/fragments.dat"
    file = open(file, "a")

    lines = all_samples(path)
    
    for line in lines:
        print(line, file=file) 