import os
import argparse
from data import all_samples

# get all observations from given directory

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-p", "--path", type=str, help="path to the directory")
    args = parser.parse_args()

    file = open(f"{os.path.dirname(args.path)}/fragments.dat", "a")

    lines = all_samples(args.path)
    for line in lines:
        print(line, file=file) 