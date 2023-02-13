import argparse
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

    aa = args.aa 
    ss = args.ss 
    range = args.range
    model = args.model 

    # read files with trained parameters
    decoder = f"{model}/decoder.h5"
    weights = f"{model}/weights.h5"
    latent = f"{model}/latent.npy"

    label = LabelMLP(aa=aa, ss=aa, r1n=range)

    # load model and predict values
    loader = DecodeLoader(decoder=decoder, weights=weights, latent=latent)
    vector = loader.predict(label.format(), 3.8)

    # change results to PDB format
    output = Output(vector=vector, bond_length=3.8)
    lines = output.to_pdb()
    for line in lines:
        print(line)