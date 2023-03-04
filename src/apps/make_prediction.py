import argparse
import logging
import numpy as np
from core.proteins import CarbonAlpha, FileParser
from core.model import DecoderLoader
from core.features import LabelMLP
from utils.structural import Vec3, Output, build_fragment, two_atoms_vector

logging.getLogger("tensorflow").disabled=True
logging.getLogger("h5py._conv").disabled=True

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-f", "--file", type=str)
    parser.add_argument("-s", "--start", type=int)
    parser.add_argument("-e", "--end", type=int)
    parser.add_argument("-aa", type=str, help="amino acids sequence")
    parser.add_argument("-ss", type=str, help="secondary structure")
    parser.add_argument("-m", "--model", type=str, help="model to be used")
    parser.add_argument("-p", "--population", type=int, help="number of fragments generated to choose the best one")
    args = parser.parse_args()

    file = args.file
    start = args.start 
    end = args.end
    model = args.model 
    population = args.population

    structure = FileParser(file=file).load_structure()

    if args.aa:
        aa = args.aa
    else:
        aa = structure.sequence(start, end)

    if args.ss:
        ss = args.ss 
    else:
        ss = ""
        for _ in range(len(aa)):
            ss += "C"

    displacement = structure.displacement(start, end)
    
    dx = displacement.x
    dy = displacement.y
    dz = displacement.z

    label = LabelMLP(aa=aa, ss=ss, dx=dx, dy=dy, dz=dz)

    decoder = DecoderLoader(decoder=f"{model}/decoder.pb", latent=f"{model}/latent.npy")

    outputs = [Output(vector=decoder.predict(label.format())[0]) for _ in range(population)]

    c_1 = structure.atoms[start-4].coordinates
    c_2 = structure.atoms[start-3].coordinates
    c_3 = structure.atoms[start-2].coordinates

    fragments = [build_fragment(c_1, c_2, c_3, output, 3.8) for output in outputs]

    lengths = [two_atoms_vector(fragment[-1], structure.atoms[end].coordinates).length() for fragment in fragments]

    errors = [np.abs(3.8 - length) for length in lengths]

    matching_fragment = fragments[errors.index(np.min(errors))] 

    for i, atom in enumerate(structure.atoms):
        index = i + 1
        if index >= start and index <= end:
            id = atom.id
            residue = atom.residue
            chain_name = atom.chain_name
            residue_id = atom.residue_id
            x = matching_fragment[index-start+3].x
            y = matching_fragment[index-start+3].y
            z = matching_fragment[index-start+3].z
            coordinates = Vec3(x=x, y=y, z=z)
            print(CarbonAlpha(id=id, residue=residue, chain_name=chain_name, residue_id=residue_id, coordinates=coordinates))
        else:
            print(atom)