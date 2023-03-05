import time
import argparse
import logging
import numpy as np
from core.model import DecoderLoader
from core.features import LabelMLP
from core.parser import FileParser, Structure, CarbonAlpha
from structural import Vec3, Output, two_atoms_vector, build_fragment

logging.getLogger("tensorflow").disabled=True
logging.getLogger("h5py._conv").disabled=True

BOND_LENGTH = 3.8

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-f", "--file", type=str, help="PDB file")
    parser.add_argument("-s", "--start", type=int, help="initial residue")
    parser.add_argument("-e", "--end", type=int, help="terminal residue")
    parser.add_argument("-m", "--model", type=str, help="model to be used")
    parser.add_argument("-p", "--population", type=int, help="number of fragments generated to choose the best one")
    args = parser.parse_args()

    file = args.file
    start = args.start 
    end = args.end 
    model = args.model 
    population = args.population

    input_structure = FileParser(file=file).load_structure() 

    aa = input_structure.read_sequence(start, end)
    ss = input_structure.read_secondary_structure(start, end)
    displacement = input_structure.local_displacement(start, end)

    dx = displacement.x
    dy = displacement.y
    dz = displacement.z

    label = LabelMLP(aa=aa, ss=ss, dx=dx, dy=dy, dz=dz)

    decoder = DecoderLoader(decoder=f"{model}/decoder.pb", latent=f"{model}/latent.npy")

    start_time = time.time()

    # raw data from decoder 
    outputs = [Output(vector=decoder.predict(label.format())[0]) for _ in range(population)]

    end_time = time.time()

    # bound atoms not included in rebuilt fragment
    c_1 = input_structure.atoms[input_structure.find_residue(start-3)].coordinates
    c_2 = input_structure.atoms[input_structure.find_residue(start-2)].coordinates
    c_3 = input_structure.atoms[input_structure.find_residue(start-1)].coordinates

    # convert generated angles to cartesian
    fragments = [build_fragment(c_1, c_2, c_3, output, BOND_LENGTH) for output in outputs] 

    lengths = []
    for fragment in fragments:
        last_bond_vector = two_atoms_vector(fragment[-1], input_structure.atoms[input_structure.find_residue(end+1)].coordinates)
        lengths.append(last_bond_vector.length())

    errors = [np.abs(BOND_LENGTH - length) for length in lengths]

    matching_fragment = fragments[errors.index(np.min(errors))] # fragment with the smallest error of bond length at last position

    new_atoms = []
    for atom in input_structure.atoms:
        atom_ss = atom.ss
        atom_id = atom.id
        atom_residue = atom.residue
        atom_chain_name = atom.chain_name
        atom_residue_id = atom.residue_id
        if atom_residue_id >= start and atom_residue_id <= end:
            vector_index = atom_residue_id - start + 3
            atom_x = matching_fragment[vector_index].x
            atom_y = matching_fragment[vector_index].y
            atom_z = matching_fragment[vector_index].z
            atom_coordinates = Vec3(x=atom_x, y=atom_y, z=atom_z)
            new_atoms.append(CarbonAlpha(ss=atom_ss, id=atom_id, residue=atom_residue, chain_name=atom_chain_name, residue_id=atom_residue_id, coordinates=atom_coordinates))
        else:
            new_atoms.append(atom)
    
    new_structure = Structure(atoms=new_atoms) # protein with inserted fragment

    lines = new_structure.to_pdb()
    for line in lines:
        print(line)

    print(f"Amino acids sequence {aa}")
    print(f"Secondary structure {ss}")
        
    total_time = end_time - start_time
    print(f"{population} outputs generated in {total_time:.3f} seconds")

    length = lengths[fragments.index(matching_fragment)]
    print(f"Last bond length {length:.3f}")