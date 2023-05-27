import os
import argparse
import logging
import numpy as np
import tensorflow as tf
from tabulate import tabulate
from core.model import DecoderLoader
from core.features import LabelMLP
from core.parser import FileParser, Structure, CarbonAlpha
from structural import Vec3, Output, two_atoms_vector, build_fragment, compute_rmsd

logging.getLogger("tensorflow").disabled=True
logging.getLogger("h5py._conv").disabled=True

BOND_LENGTH = 3.8

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-aa", type=str, help="amino acids sequence")
    parser.add_argument("-ss", type=str, help="secondary structure")
    parser.add_argument("-f", "--file", type=str, help="PDB file")
    parser.add_argument("-s", "--start", type=int, help="initial residue")
    parser.add_argument("-e", "--end", type=int, help="terminal residue")
    parser.add_argument("-m", "--model", type=str, help="model to be used")
    parser.add_argument("-r", "--repeats", type=int, help="number of returned fragments")
    parser.add_argument("-p", "--population", type=int, help="number of fragments generated to choose the best one")
    args = parser.parse_args()

    pdb = args.file
    start = args.start 
    end = args.end 
    model = args.model 
    repeats = args.repeats
    population = args.population

    input_structure = FileParser(file=pdb).load_structure() 

    if args.aa == None:
        aa = input_structure.read_sequence(start, end)
    else:
        aa = args.aa

    if args.ss == None:
        ss = input_structure.read_secondary_structure(start, end)
    else:
        ss = args.ss

    displacement = input_structure.local_displacement(start, end)

    dx = displacement.x
    dy = displacement.y
    dz = displacement.z

    decoder = DecoderLoader(decoder=f"{model}/decoder.pb", latent=f"{model}/latent.npy")
    labels = [LabelMLP(aa=aa, ss=ss, dx=dx, dy=dy, dz=dz).format() for _ in range(population)]
    labels = tf.reshape(labels, shape=(population, tf.shape(labels)[2]))

    vectors = decoder.predict(labels) # raw data from decoder
    outputs = [Output(vector) for vector in vectors]

    # bound atoms not included in rebuilt fragment
    c_1 = input_structure.atoms[input_structure.find_residue(start-3)].coordinates
    c_2 = input_structure.atoms[input_structure.find_residue(start-2)].coordinates
    c_3 = input_structure.atoms[input_structure.find_residue(start-1)].coordinates

    # convert generated angles to cartesian
    fragments = [build_fragment(c_1, c_2, c_3, output, BOND_LENGTH) for output in outputs] 

    new_structures = [] # all structures obtained from generated results
    for fragment in fragments:
        new_atoms = []
        for atom in input_structure.atoms:
            if atom.residue_id >= start and atom.residue_id <= end:
                vec_index = atom.residue_id - start + 3
                x = fragment[vec_index].x
                y = fragment[vec_index].y
                z = fragment[vec_index].z
                coordinates = Vec3(x=x, y=y, z=z)
                new_atom = CarbonAlpha(ss=atom.ss, id=atom.id, residue=atom.residue, chain_name=atom.chain_name, residue_id=atom.residue_id, coordinates=coordinates)
                new_atoms.append(new_atom)
            else:
                new_atoms.append(atom)
        
        structure = Structure(atoms=new_atoms)
        new_structures.append(structure)

    valid_structures = [] # choose structures which are not crossed 
    for structure in new_structures:
        if structure.check_if_crossing(tolerance=1.0)[0] == False:
            valid_structures.append(structure)

    last_bond_lengths = []
    for structure in valid_structures:
        last_bond_length = structure.local_displacement(end, end+1).length()
        last_bond_lengths.append(last_bond_length)

    last_bond_errors = [np.abs(BOND_LENGTH - length) for length in last_bond_lengths]
    sorted_last_bond_errors = np.sort(last_bond_errors)

    # select structures with the smallest last bond error
    if len(valid_structures) >= repeats:
        matching_structures = []
        for error in sorted_last_bond_errors[0:repeats]:
            structure = valid_structures[last_bond_errors.index(error)]
            matching_structures.append(structure)
    else:
        matching_structures = valid_structures

    pdb_name = os.path.splitext(os.path.basename(pdb))[0]
    output_path = f"{os.path.dirname(__file__)}/{pdb_name}_output.pdb"
    
    output_file = open(output_path, "a")
    
    for i, structure in enumerate(matching_structures):
        print(f"MODEL {i+1}", file=output_file)
        last_bond_vector = two_atoms_vector(structure.atoms[structure.find_residue(end)].coordinates, input_structure.atoms[input_structure.find_residue(end+1)].coordinates)
        rmsd = compute_rmsd(structure.coordinates(), input_structure.coordinates())
        print(f"{last_bond_vector.length():.3f}", file=output_file)
        print(f"{rmsd:.3f}", file=output_file)

        lines = structure.to_pdb()
        for line in lines:
            print(line, file=output_file)

    output_file.close()

    table = [["Amino acids sequence", f"{aa}"], ["Secondary structure", f"{ss}"]]
    print(tabulate(table))