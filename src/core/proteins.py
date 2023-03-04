from typing import List
from utils.structural import Vec3, two_atoms_vector

# tools for reading PDB files
# functionalities are dedicated to parse alpha carbon trace including secondary structure

RESIDUES = {"ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C", "GLN": "Q", "GLU": "E", "GLY": "G", "HIS": "H", "ILE": "I", "LEU": "L", "LYS": "K", "MET": "M", "PHE": "F", "PRO": "P", "SER": "S", "THR": "T", "TRP": "W", "TYR": "Y", "VAL": "V"}


class CarbonAlpha:
    def __init__(self, ss, id, residue, residue_id, chain_name, coordinates: Vec3):
        self._ss = ss 
        self._id = id 
        self._residue = residue
        self._residue_id = residue_id
        self._chain_name = chain_name
        self._coordinates = coordinates

    def __str__(self):
        return "%4s %6d %3s %4s %s %3d %11.3f %7.3f %7.3f %5.2f %5.2f %11s" % ("ATOM", self.id(), "CA", self.residue(), self.chain_name(), self.residue_id(), self.x(), self.y(), self.z(), 1.00, 0.00, "C")

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
    def residue_id(self):
        return self._residue_id
    
    @property
    def chain_name(self):
        return self._chain_name
    
    @property
    def coordinates(self):
        return self._coordinates
    
    @property
    def x(self):
        return self.coordinates().x 

    @property
    def y(self):
        return self.coordinates().y
    
    @property
    def z(self):
        return self.coordinates().z
    
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
    

class LineParser:
    def __init__(self, line):
        self._line = line 

    @property
    def line(self):
        return self._line 
    
    def read_id(self):
        return int(self.line()[6:11].strip())
    
    def read_residue(self):
        return self.line()[17:20]
    
    def read_residue_id(self):
        return int(self.line()[22:26].strip())
    
    def read_chain_name(self):
        return self.line()[21]
    
    def read_x(self):
        return float(self.line()[30:38].strip())
    
    def read_y(self):
        return float(self.line()[38:47].strip())
    
    def read_z(self):
        return float(self.line()[47:54].strip())
    

class FileParser:
    def __init__(self, file):
        self.file = file 
        stream = open(file)
        self._lines = stream.readlines()
        stream.close()

    @property
    def lines(self):
        return [line for line in self._lines if len(line) >= 4]
    
    def parse_ca(self):
        records = []
        for line in self.lines():
            if line[0:4] == "ATOM":
                if line[12:16].strip() == "CA":
                    records.append(line)
        return records
    
    def parse_helix(self):
        numbers = []
        for line in self.lines():
            if line[0:5] == "HELIX":
                initial = int(line[21:25].strip())
                terminal = int(line[33:37].strip())
                for i in range(initial, terminal+1):
                    numbers.append(i)
        return numbers

    def parse_sheet(self):
        numbers = []
        for line in self.lines():
            if line[0:5] == "SHEET":
                initial = int(line[22:26].strip())
                terminal = int(line[33:37].strip())
                for i in range(initial, terminal+1):
                    numbers.append(i)
        return numbers
    
    def load_atoms(self):
        atoms = []
        records = self.parse_ca()
        for record in records:
            parser = LineParser(record)
            id = parser.read_id()
            residue = parser.read_residue()
            residue_id = parser.read_residue_id()
            chain_name = parser.read_chain_name()
            x = parser.read_x()
            y = parser.read_y()
            z = parser.read_z()
            coordinates = Vec3(x=x, y=y, z=z)
            if id in self.parse_helix():
                atoms.append(CarbonAlpha(ss="H", id=id, residue=residue, residue_id=residue_id, chain_name=chain_name, coordinates=coordinates))
            if id in self.parse_sheet():
                atoms.append(CarbonAlpha(ss="E", id=id, residue=residue, residue_id=residue_id, chain_name=chain_name, coordinates=coordinates))
            else:
                atoms.append(CarbonAlpha(ss="C", id=id, residue=residue, residue_id=residue_id, chain_name=chain_name, coordinates=coordinates))
        return atoms