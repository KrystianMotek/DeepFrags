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
        # to avoid exploding gradients
        return self.alpha / 180 

    def sin_theta(self):
        return tf.sin(self.theta * np.pi / 180)

    def cos_theta(self):
        return tf.cos(self.theta * np.pi / 180)


class InputMLP(Input):
    def format(self):
        return tf.concat([self.normaliza_alpha(), self.sin_theta(), self.cos_theta()], axis=1)


class Label(ABC):
    def __init__(self, aa, ss, r1n):
        self.aa = aa
        self.ss = ss
        # distance between bound atoms
        self.r1n = r1n

    @abstractmethod
    def format(self):
        pass

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
    

class LabelMLP(Label):
    def format(self):
        return tf.concat([tf.reshape(self.r1n, shape=(1, 1)), self.encode_aa(), self.encode_ss()], axis=1)


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

    def read_alpha(self):
        ORDINAL = 9
        return tf.constant([[float(angle) for i, angle in enumerate(self.line[ORDINAL:]) if i % 2 == 0]])

    def read_theta(self):
        ORDINAL = 10
        return tf.constant([[float(angle) for i, angle in enumerate(self.line[ORDINAL:]) if i % 2 == 0]])

    def compute_r1n(self):
        ORDINAL_X = 6
        ORDINAL_Y = 7
        ORDINAL_Z = 8
        dx = float(self.line[ORDINAL_X])
        dy = float(self.line[ORDINAL_Y])
        dz = float(self.line[ORDINAL_Z])
        return tf.norm([dx, dy, dz])


class ObservationMLP(Observation):
    def create_input(self):
        return InputMLP(alpha=self.read_alpha(), theta=self.read_theta())

    def create_label(self):
        return LabelMLP(aa=self.read_aa(), ss=self.read_ss(), r1n=self.compute_r1n())


class DataSet(ABC):
    def __init__(self, file):
        self.file = file
        stream = open(file, "r")
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