import os
import warnings

warnings.filterwarnings("ignore")
os.environ['TF_CPP_MIN_LOG_LEVEL'] = '3'

import tensorflow as tf
import numpy as np
from typing import List
from abc import ABC, abstractmethod

'''
    Functionalities useful in data features extraction
    Example line containing also additional information is shown below
    4zvc A LEU 201 YKILLDSDLGS HHHHHHHCCCC -6.846 -4.808 -0.648 92.740 44.589 94.971 50.121 104.252 31.637 108.029 14.217 89.292 -127.539
'''


class Input(ABC):
    def __init__(self, alpha, theta):
        self.alpha = alpha
        self.theta = theta

    @abstractmethod
    def format(self):
        pass
    
    def normaliza_alpha(self):
        return self.alpha / 180 # to avoid exploding gradients

    def sin_theta(self):
        return tf.sin(self.theta * np.pi / 180)

    def cos_theta(self):
        return tf.cos(self.theta * np.pi / 180)


class InputMLP(Input):
    def format(self):
        return tf.concat([self.normaliza_alpha(), self.sin_theta(), self.cos_theta()], axis=1)


class Label(ABC):
    def __init__(self, aa, ss, dx, dy, dz):
        self.aa = aa
        self.ss = ss
        self.dx = dx
        self.dy = dy
        self.dz = dz

    @abstractmethod
    def format(self):
        pass

    @staticmethod
    def string_to_one_hot(string, codes):
        indices = {code: i for i, code in enumerate(codes)}
        numerical_values = [indices[value] for value in list(string)]
        depth = len(codes)
        return tf.reshape(tf.one_hot(numerical_values, depth=depth), shape=(1, depth * len(numerical_values)))
    
    @staticmethod
    def one_hot_to_string(vector, codes):
        string = ""
        categories = len(codes)
        for i in range(len(vector)):
            index = list(vector[i:i+categories]).index(1.0)
            string += codes[index]
        return string

    def encode_aa(self):
        string = self.aa
        return self.string_to_one_hot(string, codes="ARNDCQEGHILKMFPSTWYV")

    def encode_ss(self):
        string = self.ss
        return self.string_to_one_hot(string, codes="HEC")

    def displacement(self):
        dx = self.dx
        dy = self.dy
        dz = self.dz
        return tf.constant([[dx, dy, dz]])
    
    def extract_aa(self):
        pass
    
    def extract_ss(self):
        pass
    

class LabelMLP(Label):
    def format(self):
        return tf.concat([self.displacement(), self.encode_aa(), self.encode_ss()], axis=1)
    
    @staticmethod
    def extract_aa(vector):
        n = (len(vector) - 3) / 23
        ORDINAL = 20 * n + 3
        return Label.one_hot_to_string(vector[3:ORDINAL], codes="ARNDCQEGHILKMFPSTWYV")
    
    @staticmethod
    def extract_ss(vector):
        n = (len(vector) - 3) / 23
        ORDINAL = 20 * n + 3
        return Label.one_hot_to_string(vector[ORDINAL:], codes="HEC")


class Observation(ABC):
    def __init__(self, line: str):
        self.line = line.split()

    @abstractmethod
    def create_input(self) -> Input:
        pass

    @abstractmethod
    def create_label(self) -> Label:
        pass

    def read_aa(self):
        ORDINAL = 4
        length = len(self.line[ORDINAL])
        return "".join(self.line[ORDINAL][3:length-3])
    
    def read_ss(self):
        ORDINAL = 5
        length = len(self.line[ORDINAL])
        return "".join(self.line[ORDINAL][3:length-3])

    def read_dx(self):
        ORDINAL = 6
        return float(self.line[ORDINAL])

    def read_dy(self):
        ORDINAL = 7
        return float(self.line[ORDINAL])

    def read_dz(self):
        ORDINAL = 8
        return float(self.line[ORDINAL])

    def read_alpha(self):
        ORDINAL = 9
        return tf.constant([[float(angle) for i, angle in enumerate(self.line[ORDINAL:]) if i % 2 == 0]])

    def read_theta(self):
        ORDINAL = 10
        return tf.constant([[float(angle) for i, angle in enumerate(self.line[ORDINAL:]) if i % 2 == 0]])


class ObservationMLP(Observation):
    def create_input(self):
        return InputMLP(alpha=self.read_alpha(), theta=self.read_theta())

    def create_label(self):
        return LabelMLP(aa=self.read_aa(), ss=self.read_ss(), dx=self.read_dx(), dy=self.read_dy(), dz=self.read_dz())


class DataSet(ABC):
    def __init__(self, file):
        self.file = file
        stream = open(file)
        self.lines = [line for line in stream.readlines()] # read line by line
        stream.close()

    @abstractmethod
    def load_observations(self) -> List[Observation]:
        pass

    @abstractmethod
    def inputs_tensor(self):
        pass

    @abstractmethod
    def labels_tensor(self):
        pass

    def load_inputs(self):
        return [observation.create_input() for observation in self.load_observations()]

    def load_labels(self):
        return [observation.create_label() for observation in self.load_observations()]

    def save_inputs(self, file):
        np.save(file, self.inputs_tensor())

    def save_labels(self, file):
        np.save(file, self.labels_tensor())


class DataSetMLP(DataSet):
    def load_observations(self) -> List[Observation]:
        return [ObservationMLP(line) for line in self.lines]

    def inputs_tensor(self):
        return tf.concat([input.format() for input in self.load_inputs()], axis=0)

    def labels_tensor(self):
        return tf.concat([label.format() for label in self.load_labels()], axis=0)