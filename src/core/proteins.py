from typing import List
from utils.structural import Vec3, two_atoms_vector

RESIDUES = {"ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C", "GLN": "Q", "GLU": "E", "GLY": "G", "HIS": "H", "ILE": "I", "LEU": "L", "LYS": "K", "MET": "M", "PHE": "F", "PRO": "P", "SER": "S", "THR": "T", "TRP": "W", "TYR": "Y", "VAL": "V"}


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
        return "ATOM" + f"{self.id}".rjust(7) + "  " + f"CA  {self.residue} {self.chain_name}" + f"{self.residue_id}".rjust(4) + f"{self.x:.3f}".rjust(12) + f"{self.y:.3f}".rjust(8) + f"{self.z:.3f}".rjust(8)
    
    def translate_three_letters(self):
        return RESIDUES[self.residue]


class Record:
    def __init__(self, line: str):
        self.line = line
    
    def read_id(self):
        ORDINAL = 1
        return self.line.split()[ORDINAL]
    
    def read_residue(self):
        ORDINAL = 3
        return self.line.split()[ORDINAL]
    
    def read_chain_name(self):
        ORDINAL = 4
        return self.line.split()[ORDINAL]
    
    def read_residue_id(self):
        ORDINAL = 5
        return self.line.split()[ORDINAL]
    
    def read_coordinates(self):
        ORDINAL_X = 6
        ORDINAL_Y = 7
        ORDINAL_Z = 8
        x = float(self.line.split()[ORDINAL_X])
        y = float(self.line.split()[ORDINAL_Y])
        z = float(self.line.split()[ORDINAL_Z])
        return Vec3(x=x, y=y, z=z)
    
    def transform(self):
        return CarbonAlpha(id=self.read_id(), residue=self.read_residue(), chain_name=self.read_chain_name(), residue_id=self.read_residue_id(), coordinates=self.read_coordinates())


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

    def select_lines(self) -> List[Record]:
        records = []
        for line in self.lines:
            if len(line.strip()) != 0:
                if line.split()[0] == "ATOM" and line.split()[2] == "CA":
                    records.append(Record(line=line))
        return records
    
    def load_structure(self):
        return Structure(atoms=[record.transform() for record in self.select_lines()])