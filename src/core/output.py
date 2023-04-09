import math 
import tensorflow as tf 


class Output:
    def __init__(self, vector):
        self._vector = vector 

    @property
    def vector(self):
        return self._vector
    
    def n(self):
        # number of amino acids
        return int(len(self.vector) / 3)

    def alpha(self):
        start = 0
        end = self.n()
        return [self.vector[i] * 180.0 for i in range(start, end)]

    def sin_theta(self):
        start = self.n()
        end = 2 * self.n()
        return [self.vector[i] for i in range(start, end)]

    def cos_theta(self):
        start = 2 * self.n()
        end = 3 * self.n()
        return [self.vector[i] for i in range(start, end)]

    def theta(self):
        sin_theta = self.sin_theta()
        cos_theta = self.cos_theta()
        return [sin_cos_to_angle(sin, cos) for sin, cos in zip(sin_theta, cos_theta)]
    

def to_degrees(radians):
    return radians * 180.0 / math.pi 


def to_radians(degrees):
    return degrees * math.pi / 180.0


def sin_cos_to_angle(sin, cos):
    k = math.sqrt(sin ** 2 + cos ** 2)
    return to_degrees(math.atan2(sin / k, cos / k))


def angles_to_cartesian(atom_1, atom_2, atom_3, bond_length, alpha, theta):
    sin_alpha = tf.sin(to_radians(alpha))
    cos_alpha = tf.cos(to_radians(alpha))
    sin_theta = tf.sin(to_radians(theta))
    cos_theta = tf.cos(to_radians(theta))

    x = bond_length * cos_alpha
    y = bond_length * sin_alpha * cos_theta
    z = bond_length * sin_alpha * sin_theta

    v_12 = tf.subtract(atom_2, atom_1)
    v_23 = tf.subtract(atom_3, atom_2)
    
    k = tf.linalg.cross(v_12, tf.nn.l2_normalize(v_23))

    l = tf.linalg.cross(tf.nn.l2_normalize(k), tf.nn.l2_normalize(v_23))

    new_x = atom_3[0] - tf.nn.l2_normalize(v_23)[0] * x + l[0] * y + tf.nn.l2_normalize(k)[0] * z 
    new_y = atom_3[1] - tf.nn.l2_normalize(v_23)[1] * x + l[1] * y + tf.nn.l2_normalize(k)[1] * z 
    new_z = atom_3[2] - tf.nn.l2_normalize(v_23)[2] * x + l[2] * y + tf.nn.l2_normalize(k)[2] * z 

    return tf.constant([new_x.numpy(), new_y.numpy(), new_z.numpy()])


def build_fragment(c_1, c_2, c_3, output: Output, bond_length):
    n = output.n()
    alpha = output.alpha()
    theta = output.theta()

    start = 3
    end = n + 3
    atoms = [c_1, c_2, c_3]
    for i in range(start, end):
        c_i = atoms[i-3]
        c_j = atoms[i-2]
        c_k = atoms[i-1]
        c_new = angles_to_cartesian(c_i, c_j, c_k, bond_length, alpha[i-3], theta[i-3])
        atoms.append(c_new)

    return atoms[3:]