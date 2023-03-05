from typing import List
from structural import Vec3, two_atoms_vector

# tools for reading PDB files
# functionalities are dedicated to parse alpha carbon trace including secondary structure

RESIDUES = {"ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C", "GLN": "Q", "GLU": "E", "GLY": "G", "HIS": "H", "ILE": "I",
            "LEU": "L", "LYS": "K", "MET": "M", "PHE": "F", "PRO": "P", "SER": "S", "THR": "T", "TRP": "W", "TYR": "Y", "VAL": "V"}
 

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

    def find_atom(self, id):
        for atom in self.atoms:
            if atom.residue_id == id:
                return self.atoms.index(atom)
            else:
                continue

    def sequence(self, i, j):
        string = ""
        atoms = self.atoms[self.find_atom(i):self.find_atom(j)+1]
        for atom in atoms:
            string += RESIDUES[atom.residue]
        return string 
    
    def secondary_structure(self, i, j):
        string = ""
        atoms = self.atoms[self.find_atom(i):self.find_atom(j)+1]
        for atom in atoms:
            string += atom.ss
        return string
    
    def displacement(self, i, j):
        coordinates_i = self.atoms[self.find_atom(i)].coordinates
        coordinates_j = self.atoms[self.find_atom(j)].coordinates
        return two_atoms_vector(coordinates_i, coordinates_j)
    
    def to_pdb(self):
        return [atom.__str__() for atom in self.atoms]
    

class LineParser:
    def __init__(self, line):
        self._line = line 

    @property
    def line(self):
        return self._line 
    
    def read_id(self):
        START = 6
        END = 11
        return int(self.line[START:END].strip())
    
    def read_residue(self):
        START = 17
        END = 20
        return self.line[START:END]
    
    def read_residue_id(self):
        START = 22
        END = 26
        return int(self.line[START:END].strip())
    
    def read_chain_name(self):
        ORDINAL = 21
        return self.line[ORDINAL]
    
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
            START_ATOM = 0
            END_ATOM = 4
            if line[START_ATOM:END_ATOM] == "ATOM":
                START_CA = 12
                END_CA = 16
                if line[START_CA:END_CA].strip() == "CA":
                    records.append(line)
        return records
    
    def parse_helix(self):
        numbers = []
        for line in self.lines:
            START_HELIX = 0
            END_HELIX = 5
            if line[START_HELIX:END_HELIX] == "HELIX":
                START_INITIAL = 21
                END_INITIAL = 25
                initial = int(line[START_INITIAL:END_INITIAL].strip())
                START_TERMINAL = 33
                END_TERMINAL = 37
                terminal = int(line[START_TERMINAL:END_TERMINAL].strip())
                for i in range(initial, terminal+1):
                    numbers.append(i)
        return numbers

    def parse_sheet(self):
        numbers = []
        for line in self.lines:
            START_SHEET = 0
            END_SHEET = 5
            if line[START_SHEET:END_SHEET] == "SHEET":
                START_INITIAL = 22
                END_INITIAL = 26
                initial = int(line[START_INITIAL:END_INITIAL].strip())
                START_TERMINAL = 33
                END_TERMINAL = 37
                terminal = int(line[START_TERMINAL:END_TERMINAL].strip())
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
            if id not in self.parse_helix() and id not in self.parse_sheet():
                atoms.append(CarbonAlpha(ss="C", id=id, residue=residue, residue_id=residue_id, chain_name=chain_name, coordinates=coordinates))
            else:
                if id in self.parse_helix():
                    atoms.append(CarbonAlpha(ss="H", id=id, residue=residue, residue_id=residue_id, chain_name=chain_name, coordinates=coordinates))
                if id in self.parse_sheet():
                    atoms.append(CarbonAlpha(ss="E", id=id, residue=residue, residue_id=residue_id, chain_name=chain_name, coordinates=coordinates))
        return atoms

    def load_structure(self):
        return Structure(atoms=self.load_atoms())