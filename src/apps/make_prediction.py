import argparse
import logging
import numpy as np
from core.proteins import CarbonAlpha, Structure, FileParser
from core.model import DecoderLoader
from core.features import LabelMLP
from structural import Vec3, Output, build_fragment, two_atoms_vector

logging.getLogger("tensorflow").disabled=True
logging.getLogger("h5py._conv").disabled=True

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-f", "--file", type=str)
    parser.add_argument("-s", "--start", type=int)
    parser.add_argument("-e", "--end", type=int)
    parser.add_argument("-m", "--model", type=str, help="model to be used")
    parser.add_argument("-p", "--population", type=int, help="number of fragments generated to choose the best one")
    args = parser.parse_args()

    file = args.file
    start = args.start 
    end = args.end
    model = args.model 
    population = args.population

    structure = FileParser(file=file).load_structure() 

    aa = structure.sequence(start, end)
    ss = structure.secondary_structure(start, end)

    displacement = structure.displacement(start, end) 

    dx = displacement.x
    dy = displacement.y
    dz = displacement.z

    label = LabelMLP(aa=aa, ss=ss, dx=dx, dy=dy, dz=dz)

    decoder = DecoderLoader(decoder=f"{model}/decoder.pb", latent=f"{model}/latent.npy")

    # raw data from decoder 
    outputs = [Output(vector=decoder.predict(label.format())[0]) for _ in range(population)]

    # boundary atoms not included in rebuilt fragment
    c_1 = structure.atoms[start-4].coordinates
    c_2 = structure.atoms[start-3].coordinates
    c_3 = structure.atoms[start-2].coordinates

    # convert generated angles to cartesian
    fragments = [build_fragment(c_1, c_2, c_3, output, 3.8) for output in outputs] 

    lengths = [two_atoms_vector(fragment[-1], structure.atoms[end].coordinates).length() for fragment in fragments] 

    errors = [np.abs(3.8 - length) for length in lengths]

    matching_fragment = fragments[errors.index(np.min(errors))] # fragment with the smallest error of bond length at last position

    new_atoms = []
    for i, atom in enumerate(structure.atoms):
        if i >= start and i <= end:
            ss = atom.ss
            id = atom.id
            residue = atom.residue
            chain_name = atom.chain_name
            residue_id = atom.residue_id
            x = matching_fragment[i-start+3].x
            y = matching_fragment[i-start+3].y
            z = matching_fragment[i-start+3].z
            coordinates = Vec3(x=x, y=y, z=z)
            new_atoms.append(CarbonAlpha(ss, id=id, residue=residue, chain_name=chain_name, residue_id=residue_id, coordinates=coordinates))
        else:
            new_atoms.append(atom)

    # protein with inserted fragment
    new_structure = Structure(atoms=new_atoms)

    lines = new_structure.to_pdb()
    for line in lines:
        print(line)