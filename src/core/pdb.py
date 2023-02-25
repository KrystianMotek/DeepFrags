from typing import List
from structural import Vec3, two_atoms_vector

RESIDUES = {"A": "ALA", "R": "ARG", "N": "ASN", "D": "ASP", "C": "CYS", "Q": "GLN", "E": "GLU", "G": "GLY", "H": "HIS", "I": "ILE", "L": "LEU", "K": "LYS", "M": "MET", "F": "PHE", "P": "PRO", "S": "SER", "T": "THR", "W": "TRP", "Y": "TYR", "V": "VAL"}


class Atom:
    def __init__(self, id, name, residue, chain_name, residue_id, coordinates: Vec3):
        self.id = id
        self.name = name
        self.residue = residue
        self.chain_name = chain_name
        self.residue_id = residue_id
        self.coordinates = coordinates

        self.x = self.coordinates.x
        self.y = self.coordinates.y
        self.z = self.coordinates.z

    def check_if_ca(self):
        return self.name == "CA"
    
    def to_pdb(self):
        pass


class Line:
    def __init__(self, line: str):
        self.line = line.split()

    def read_id(self):
        ORDINAL = 1
        return self.line[ORDINAL]

    def read_name(self):
        ORDINAL = 2
        return self.line[ORDINAL]
    
    def read_residue(self):
        ORDINAL = 3
        return self.line[ORDINAL]
    
    def read_chain_name(self):
        ORDINAL = 4
        return self.line[ORDINAL]
    
    def read_residue_id(self):
        ORDINAL = 5
        return self.line[ORDINAL]

    def read_coordinates(self):
        ORDINAL_X = 6
        ORDINAL_Y = 7
        ORDINAL_Z = 8
        x = float(self.line[ORDINAL_X])
        y = float(self.line[ORDINAL_Y])
        z = float(self.line[ORDINAL_Z])
        return Vec3(x=x, y=y, z=z)
    
    def load_atom(self):
        return Atom(id=self.read_id(), name=self.read_name(), residue=self.read_residue(), chain_name=self.read_chain_name(), residue_id=self.read_residue_id(), coordinates=self.read_coordinates())
    

class Structure:
    def __init__(self, atoms: List[Atom]):
        self.atoms = atoms

    def extract_ca(self):
        self.atoms = [atom for atom in self.atoms if atom.check_if_ca()]
    
    def local_distance(self, index_1, index_2):
        atom_1 = self.atoms[index_1]
        atom_2 = self.atoms[index_2]
        vector = two_atoms_vector(atom_1, atom_2)
        return vector.length()

    def to_pdb(self):
        pass


class File:
    def __init__(self, path):
        self.stream = open(path, "r")
        self.lines = self.stream.readlines()