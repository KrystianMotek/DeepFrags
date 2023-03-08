import argparse
import numpy as np
from core.model import DecoderLoader
from core.features import LabelMLP
from structural import Output

BOND_LENGTH = 3.8

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-l", "--labels", type=str)
    parser.add_argument("-m", "--model", type=str)
    args = parser.parse_args()
    
    labels = np.load(args.labels)

    ss = []
    for label in labels:
        ss.append(LabelMLP.extract_ss(label))

    model = args.model 
    decoder = DecoderLoader(decoder=f"{model}/decoder.pb", latent=f"{model}/latent.npy")

    reconstructed_data = decoder.predict(labels)

    alpha = []
    theta = []
    for vector in reconstructed_data:
        output = Output(vector=vector)
        alpha.append(output.alpha())
        theta.append(output.theta())

    # planar and dihedral angles correlation
    # histograms for each secondary structure