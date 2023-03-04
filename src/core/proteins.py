from typing import List
from utils.structural import Vec3, two_atoms_vector

# tools for reading PDB files

RESIDUES = {"ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C", "GLN": "Q", "GLU": "E", "GLY": "G", "HIS": "H", "ILE": "I",
            "LEU": "L", "LYS": "K", "MET": "M", "PHE": "F", "PRO": "P", "SER": "S", "THR": "T", "TRP": "W", "TYR": "Y", "VAL": "V"}


class CarbonAlpha:
    def __init__(self, ss, id, residue, chain_name, coordinates: Vec3):
        self._ss = ss 
        self._id = id 
        self._residue = residue
        self._chain_name = chain_name
        self._coordinates = coordinates

    def __str__(self):
        pass

    @property
    def ss(self):
        return self._ss 

    @property
    def id(self):
        return self._id 

    @property
    def residue(self):
        return self._residue
    
    @property
    def chain_name(self):
        return self._chain_name
    
    @property
    def coordinates(self):
        return self._coordinates
    
    @coordinates.setter
    def coordinates(self, vector: Vec3):
        self._coordinates = vector


class Chain:
    def __init__(self, atoms: List[CarbonAlpha]):
        self._atoms = atoms 
    
    @property
    def atoms(self):
        return self._atoms

    def sequence(self, i, j):
        string = ""
        atoms = self.atoms()[i-1:j-1]
        for atom in atoms:
            string += RESIDUES[atom.residue()]
        return string 
    
    def secondary_structure(self, i, j):
        string = ""
        atoms = self.atoms()[i-1:j-1]
        for atom in atoms:
            string += atom.ss()
        return string
    
    def displacement(self, i, j):
        coordinates_i = self.atoms[i-1].coordinates()
        coordinates_j = self.atoms[j-1].coordinates()
        vector = two_atoms_vector(coordinates_i, coordinates_j)
        return vector.length()


class Structure:
    def __init__(self, chains: List[Chain]):
        self._chains = chains

    @property
    def chains(self):
        return self._chains
    
    def to_pdb(self):
        lines = []
        for chain in self.chains():
            for atom in chain.atoms():
                lines.append(atom.__str__())
        return lines