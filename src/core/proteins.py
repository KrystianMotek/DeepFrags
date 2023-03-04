from typing import List
from utils.structural import Vec3, two_atoms_vector

# tools for reading PDB files

RESIDUES = {"ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C", "GLN": "Q", "GLU": "E", "GLY": "G", "HIS": "H", "ILE": "I", "LEU": "L", "LYS": "K", "MET": "M", "PHE": "F", "PRO": "P", "SER": "S", "THR": "T", "TRP": "W", "TYR": "Y", "VAL": "V"}


class Atom:
    def __init__(self, id, name, coordinates: Vec3):
        self.id = id 
        self.name = name 
        self.coordinates = coordinates

        self.x = self.coordinates.x
        self.y = self.coordinates.y
        self.z = self.coordinates.z

    def element(self):
        if self.name in ["NA", "K", "MG", "ZN", "CU", "FE", "SE", "CO", "MN", "MO"]:
            return self.name 
        else:
            return self.name[0]

    def is_ca(self):
        return self.name == "CA"


class Residue:
    def __init__(self, id, name, atoms: List[Atom]):
        self.id = id 
        self.name = name
        self.atoms = atoms
    
    def translate_three_letters(self):
        return RESIDUES[self.name]


class Chain:
    def __init__(self, name, residues: List[Residue]):
        self.name = name
        self.residues = residues

    def sequence(self, i, j):
        string = ""
        for residue in self.residues[i-1:j-1]:
            string += residue.translate_three_letters()
        return string
    
    def extract_ca(self):
        trace = []
        for residue in self.residues:
            for atom in residue:
                if atom.is_ca():
                    trace.append(atom)
        return trace

    def displacement(self):
        pass


class Structure:
    def __init__(self, chains: List[Chain]):
        self.chain = chains


class LineParser:
    def __init__(self, line: str):
        self.line = line
    
    def read_id(self):
        START = 6
        END = 11
        return int(self.line[START:END].strip())

    def read_name(self):
        START = 12
        END = 16
        return self.line[START:END].strip()
    
    def read_x(self):
        START = 30
        END = 38
        return float(self.line[START:END].strip())
    
    def read_y(self):
        START = 38
        END = 47
        return float(self.line[START:END].strip())
    
    def read_z(self):
        START = 47
        END = 54
        return float(self.line[START:END].strip())
    
    def load_atom(self):
        return Atom(id=self.read_id(), name=self.read_name(), coordinates=Vec3(x=self.read_x(), y=self.read_y(), z=self.read_z()))
    

class FileParser:
    def __init__(self, file):
        self.file = file