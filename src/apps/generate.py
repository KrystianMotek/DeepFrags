import argparse
import random
from data import all_samples

# random observations selection from all the files in given directory

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-p", "--path", type=str, help="path to the directory")
    parser.add_argument("-s", "--size", type=int, help="number of samples")
    args = parser.parse_args()

    path = args.path
    size = args.size
    work_file = open(f"{path}/fragments.dat", "a")

    # read each file and choose samples
    lines = all_samples(args.path)
    lines = random.sample(lines, size)

    for line in lines:
        print(line, file=work_file) 