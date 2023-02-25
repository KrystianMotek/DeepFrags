from typing import List
from structural import Vec3, two_atoms_vector

RESIDUES = {"A": "ALA", "R": "ARG", "N": "ASN", "D": "ASP", "C": "CYS", "Q": "GLN", "E": "GLU", "G": "GLY", "H": "HIS", "I": "ILE", "L": "LEU", "K": "LYS", "M": "MET", "F": "PHE", "P": "PRO", "S": "SER", "T": "THR", "W": "TRP", "Y": "TYR", "V": "VAL"}


class CarbonAlpha:
    def __init__(self, id, residue, chain_name, residue_id, coordinates):
        self.id = id
        self.residue = residue
        self.chain_name = chain_name
        self.residue_id = residue_id
        self.coordinates = coordinates

        self.x = self.coordinates.x
        self.y = self.coordinates.y
        self.z = self.coordinates.z

    def __str__(self):
        pass


class Structure:
    def __init__(self, atoms):
        self.atoms = atoms
    
    def to_pdb(self):
        return [atom.__str__() for atom in self.atoms]
    
    def local_distance(self, i, j):
        coordinates_i = self.atoms[i-1].coordinates
        coordinates_j = self.atoms[j-1].coordinates
        vector = two_atoms_vector(coordinates_i, coordinates_j)
        return vector.length()


class FileParser:
    def __init__(self, file):
        self.file = file
        self.stream = open(self.file)
        self.lines = self.stream.readlines()

    def load_structure(self):
        pass