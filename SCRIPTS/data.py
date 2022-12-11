import argparse
import warnings
import os

warnings.filterwarnings("ignore")
os.environ["TF_CPP_MIN_LOG_LEVEL"] = "3"

import tensorflow as tf
import numpy as np
from statistics.statistics import hec_distribution

'''
    Example line from data file
    1a62A 3 THR 4 MNLTELKNTPV CCHHHHHCCCH 1.284 -5.477 -2.600 99.739 -124.579 91.327 49.439 90.431 53.359 90.948 45.020 87.256 56.852
    Observation class represents a single set of features 
    DataSet class holds data from all observations for a given number of amino acids
'''


class Input:
    def __init__(self, displacement, alpha, theta):
        self.displacement = displacement
        self.alpha = alpha
        self.theta = theta

    def normalize_alpha(self):
        # to avoid exploding gradients
        return self.alpha / 180
    
    def sin_theta(self):
        return tf.sin(self.theta * np.pi / 180)
    
    def cos_theta(self):
        return tf.cos(self.theta * np.pi / 180)
    
    def to_vector(self):
        return tf.concat([self.displacement, self.normalize_alpha(), self.sin_theta(), self.cos_theta()], axis=1)

    
class Label:
    def __init__(self, aa, ss):
        self.aa = aa
        self.ss = ss

    @staticmethod
    def string_to_one_hot(string, codes):
        indices = {code: i for i, code in enumerate(codes)}
        numerical_values = [indices[value] for value in list(string)]
        depth = len(codes)
        return tf.reshape(tf.one_hot(numerical_values, depth=depth), shape=(1, depth * len(numerical_values)))

    def encode_aa(self):
        string = self.aa
        return self.string_to_one_hot(string, codes="ARNDCQEGHILKMFPSTWYV")

    def encode_ss(self):
        string = self.ss
        return self.string_to_one_hot(string, codes="HEC")
    
    def to_vector(self):
        return tf.concat([self.encode_aa(), self.encode_ss()], axis=1)
    

class Observation:
    def __init__(self, line):
        self.line = line # line should be featured as a list

    def aa(self):
        ORDINAL = 4
        length = len(self.line[ORDINAL])
        return "".join(self.line[ORDINAL][3:length-3])
    
    def ss(self):
        ORDINAL = 5
        length = len(self.line[ORDINAL])
        return "".join(self.line[ORDINAL][3:length-3])
    
    def displacement(self):
        ORDINAL_X = 6
        ORDINAL_Y = 7
        ORDINAL_Z = 8
        dx = float(self.line[ORDINAL_X])
        dy = float(self.line[ORDINAL_Y])
        dz = float(self.line[ORDINAL_Z])
        return tf.constant([[dx, dy, dz]])
    
    def alpha(self):
        ORDINAL = 9
        return tf.constant([[float(angle) for i, angle in enumerate(self.line[ORDINAL:]) if i % 2 == 0]])

    def theta(self):
        ORDINAL = 10
        return tf.constant([[float(angle) for i, angle in enumerate(self.line[ORDINAL:]) if i % 2 == 0]])
    
    def load_input(self):
        coordinate = Input(self.displacement(), self.alpha(), self.theta())
        return coordinate.to_vector()
    
    def load_label(self):
        label = Label(self.aa(), self.ss())
        return label.to_vector()
    

class DataSet:
    def __init__(self, file):
        self.file = file
        self.lines = [line.split() for line in open(file, "r").readlines()]
        self.observations = [Observation(line) for line in self.lines]

    def inputs_tensor(self):
        return tf.concat([observation.load_input() for observation in self.observations], axis=0)
    
    def labels_tensor(self):
        return tf.concat([observation.load_label() for observation in self.observations], axis=0)
    
    def save_inputs(self, file):
        np.save(file, self.inputs_tensor())

    def save_labels(self, file):
        np.save(file, self.labels_tensor())


if __name__ == "__main__":
    tf.compat.v1.logging.set_verbosity(tf.compat.v1.logging.ERROR)
    
    parser = argparse.ArgumentParser()
    parser.add_argument("-f", "--file", type=str, help="path to the data file")
    args = parser.parse_args()

    work_directory = os.path.dirname(os.path.abspath(__file__))

    # load and save
    data = DataSet(args.file)
    data.save_inputs(f"{work_directory}/inputs")
    data.save_labels(f"{work_directory}/labels")

    observations = data.observations

    # get statistics of plane and dihedral angles in data set
    ss_alpha_theta = [[observation.ss(), observation.alpha(), observation.theta()] for observation in observations]
    angles = np.concatenate([[[ss, float(alpha), float(theta)] for ss, alpha, theta in zip(element[0], element[1][0], element[2][0])] for element in ss_alpha_theta])

    # each list represents different secondary structure
    h_angles, e_angles, c_angles = hec_distribution(angles)[0], hec_distribution(angles)[1], hec_distribution(angles)[2]
