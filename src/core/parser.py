from typing import List
from structural import Vec3, two_atoms_vector

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
        formatters = ("ATOM", self.id, "CA", self.residue, self.chain_name, self.residue_id, self.x, self.y, self.z, 1.00, 0.00, "C") 
        return "%4s %6d %3s %4s %s %3d %11.3f %7.3f %7.3f %5.2f %5.2f %11s" % formatters
    
    def __lt__(self, other):
        return self.residue_id < other.residue_id 
    
    def __gt__(self, other):
        return self.residue_id > other.residue_id

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
        return self.coordinates.x 

    @property
    def y(self):
        return self.coordinates.y
    
    @property
    def z(self):
        return self.coordinates.z
    
    @coordinates.setter
    def coordinates(self, vector: Vec3):
        self._coordinates = vector


class Structure:
    def __init__(self, atoms: List[CarbonAlpha]):
        self._atoms = atoms 
    
    @property
    def atoms(self):
        return self._atoms
    
    def coordinates(self):
        # get list of coordinates of all atoms
        return [atom.coordinates for atom in self.atoms]

    def find_residue(self, residue_id):
        for atom in self.atoms:
            if atom.residue_id == residue_id:
                return self.atoms.index(atom)
            else:
                continue

    def read_sequence(self, i, j):
        string = ""
        atoms = self.atoms[self.find_residue(i):self.find_residue(j)+1]
        for atom in atoms:
            string += RESIDUES[atom.residue]
        return string 
    
    def read_secondary_structure(self, i, j):
        string = ""
        atoms = self.atoms[self.find_residue(i):self.find_residue(j)+1]
        for atom in atoms:
            string += atom.ss
        return string 
    
    def local_displacement(self, i, j):
        coordinates_i = self.atoms[self.find_residue(i)].coordinates
        coordinates_j = self.atoms[self.find_residue(j)].coordinates
        return two_atoms_vector(coordinates_i, coordinates_j)
    
    def to_pdb(self):
        return [atom.__str__() for atom in self.atoms]
    

class LineParser:
    def __init__(self, line):
        self._line = line 

    @property
    def line(self):
        return self._line 
    
    def parse_id(self):
        ORDINAL_START = 6
        ORDINAL_END = 11
        return int(self.line[ORDINAL_START:ORDINAL_END].strip())
    
    def parse_residue(self):
        ORDINAL_START = 17
        ORDINAL_END = 20
        return self.line[ORDINAL_START:ORDINAL_END]
    
    def parse_residue_id(self):
        ORDINAL_START = 22
        ORDINAL_END = 26
        return int(self.line[ORDINAL_START:ORDINAL_END].strip())
    
    def parse_chain_name(self):
        ORDINAL = 21
        return self.line[ORDINAL]
    
    def parse_x(self):
        ORDINAL_START = 30
        ORDINAL_END = 38
        return float(self.line[ORDINAL_START:ORDINAL_END].strip())
    
    def parse_y(self):
        ORDINAL_START = 38
        ORDINAL_END = 47
        return float(self.line[ORDINAL_START:ORDINAL_END].strip())
    
    def parse_z(self):
        ORDINAL_START = 47
        ORDINAL_END = 54
        return float(self.line[ORDINAL_START:ORDINAL_END].strip())
    

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
        for line in self.lines:
            ORDINAL_START_ATOM = 0
            ORDINAL_END_ATOM = 4
            if line[ORDINAL_START_ATOM:ORDINAL_END_ATOM] == "ATOM":
                ORDINAL_START_CA = 12
                ORDINAL_END_CA = 16
                if line[ORDINAL_START_CA:ORDINAL_END_CA].strip() == "CA":
                    records.append(line)
        return records
    
    def parse_helix(self):
        numbers = []
        for line in self.lines:
            ORDINAL_START_HELIX = 0
            ORDINAL_END_HELIX = 5
            if line[ORDINAL_START_HELIX:ORDINAL_END_HELIX] == "HELIX":
                ORDINAL_START_INITIAL = 21
                ORDINAL_END_INITIAL = 25
                initial_residue_id = int(line[ORDINAL_START_INITIAL:ORDINAL_END_INITIAL].strip())
                ORDINAL_START_TERMINAL = 33
                ORDINAL_END_TERMINAL = 37
                terminal_residue_id = int(line[ORDINAL_START_TERMINAL:ORDINAL_END_TERMINAL].strip())
                for i in range(initial_residue_id, terminal_residue_id+1):
                    numbers.append(i)
        return numbers

    def parse_sheet(self):
        numbers = []
        for line in self.lines:
            ORDINAL_START_SHEET = 0
            ORDINAL_END_SHEET = 5
            if line[ORDINAL_START_SHEET:ORDINAL_END_SHEET] == "SHEET":
                ORDINAL_START_INITIAL = 22
                ORDINAL_END_INITIAL = 26
                initial_residue_id = int(line[ORDINAL_START_INITIAL:ORDINAL_END_INITIAL].strip())
                ORDINAL_START_TERMINAL = 33
                ORDINAL_END_TERMINAL = 37
                terminal_residue_id = int(line[ORDINAL_START_TERMINAL:ORDINAL_END_TERMINAL].strip())
                for i in range(initial_residue_id, terminal_residue_id+1):
                    numbers.append(i)
        return numbers
    
    def load_atoms(self):
        atoms = []
        records = self.parse_ca()
        for record in records:
            parser = LineParser(record)
            id = parser.parse_id()
            residue = parser.parse_residue()
            residue_id = parser.parse_residue_id()
            chain_name = parser.parse_chain_name()
            x = parser.parse_x()
            y = parser.parse_y()
            z = parser.parse_z()
            coordinates = Vec3(x=x, y=y, z=z)
            # search secondary structure for each atom
            if residue_id not in self.parse_helix() and residue_id not in self.parse_sheet():
                atoms.append(CarbonAlpha(ss="C", id=id, residue=residue, residue_id=residue_id, chain_name=chain_name, coordinates=coordinates))
            else:
                if residue_id in self.parse_helix():
                    atoms.append(CarbonAlpha(ss="H", id=id, residue=residue, residue_id=residue_id, chain_name=chain_name, coordinates=coordinates))
                if residue_id in self.parse_sheet():
                    atoms.append(CarbonAlpha(ss="E", id=id, residue=residue, residue_id=residue_id, chain_name=chain_name, coordinates=coordinates))
        return atoms

    def load_structure(self):
        return Structure(atoms=self.load_atoms())