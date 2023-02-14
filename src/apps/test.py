import argparse
import logging
from model import DecodeLoader
from features import LabelMLP
from statistical import Output

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-aa", type=str, help="amino acids sequence")
    parser.add_argument("-ss", type=str, help="secondary structure")
    parser.add_argument("-r", "--range", type=float, help="distance from first atom to the last one")
    parser.add_argument("-m", "--model", type=str, help="model to be used")
    args = parser.parse_args()

    logging.getLogger("tensorflow").disabled=True
    logging.getLogger("h5py._conv").disabled=True

    aa = args.aa 
    ss = args.ss 
    range = args.range
    model = args.model 

    decoder = f"{model}/decoder.h5"
    latent = f"{model}/latent.npy"

    label = LabelMLP(aa=aa, ss=aa, r1n=range)

    # load model and predict values
    loader = DecodeLoader(decoder=decoder, latent=latent)
    vector = loader.predict(label.format(), 3.8)

    # change results to PDB format
    output = Output(vector=vector, bond_length=3.8)
    lines = output.to_pdb()
    for line in lines:
        print(line)