import argparse
import logging
import numpy as np
from proteins import FileParser
from model import DecoderLoader
from features import LabelMLP
from utils.structural import Output

logging.getLogger("tensorflow").disabled=True
logging.getLogger("h5py._conv").disabled=True

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-f", "--file", type=str, help="PDB file")
    parser.add_argument("-s", "--start", type=int)
    parser.add_argument("-e", "--end", type=int)
    parser.add_argument("-m", "--model", type=str, help="model to be used")
    parser.add_argument("-p", "--population", type=int, help="number of fragments generated to choose the best one")
    parser.add_argument("-i", "--inverse", default=False, action=argparse.BooleanOptionalAction, help="building from the last atom")
    args = parser.parse_args()

    file = args.file
    start = args.start 
    end = args.end
    model = args.model
    population = args.population

    decoder = f"{model}/decoder.pb"
    latent = f"{model}/latent.npy"

    if args.inverse:
        pass
    else:
        pass