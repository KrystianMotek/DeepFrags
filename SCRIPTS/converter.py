import argparse
import warnings
import os

warnings.filterwarnings("ignore")
os.environ["TF_CPP_MIN_LOG_LEVEL"] = "3"

import tensorflow as tf
import numpy as np
from data.data import DataSet
from statistics.statistics import hec_distribution
from statistics.graphs import alpha_histogram, theta_histogram, correlation_plot

'''
    This script reads data from file and converts to binary format
    Moreover statistics of given data can be automatically generating
    Obtained plots are stored in a dedicated directory
'''

if __name__ == "__main__":
    tf.compat.v1.logging.set_verbosity(tf.compat.v1.logging.ERROR)
    
    parser = argparse.ArgumentParser()
    parser.add_argument("-f", "--file", type=str, help="path to the data file")
    parser.add_argument("-g", "--graphs", type=eval, choices=[True, False], help="make plots choice", default=False)
    args = parser.parse_args()
    
    # directory where generated files will be stored
    work_directory = os.path.dirname(args.file)

    # load and save
    data = DataSet(args.file)
    data.save_inputs(f"{work_directory}/inputs")
    data.save_labels(f"{work_directory}/labels")

    if args.graphs == True:
        observations = data.observations

        # get statistics of plane and dihedral angles in data set
        ss_alpha_theta = [[observation.ss(), observation.alpha(), observation.theta()] for observation in observations]
        angles = np.concatenate([[[ss, float(alpha), float(theta)] for ss, alpha, theta in zip(element[0], element[1][0], element[2][0])] for element in ss_alpha_theta])

        # each list represents different secondary structure
        h_angles, e_angles, c_angles = hec_distribution(angles)[0], hec_distribution(angles)[1], hec_distribution(angles)[2]

        h_alpha = [float(angle[0]) for angle in h_angles]
        h_theta = [float(angle[1]) for angle in h_angles]

        e_alpha = [float(angle[0]) for angle in e_angles]
        e_theta = [float(angle[1]) for angle in e_angles]

        c_alpha = [float(angle[0]) for angle in c_angles]
        c_theta = [float(angle[1]) for angle in c_angles]

        graphs = f"{work_directory}/graphs" 
        os.makedirs(graphs)

        alpha_histogram(h_alpha, f"{graphs}/h_alpha")
        theta_histogram(h_theta, f"{graphs}/h_theta")

        alpha_histogram(e_alpha, f"{graphs}/e_alpha")
        theta_histogram(e_theta, f"{graphs}/e_theta")

        alpha_histogram(c_alpha, f"{graphs}/c_alpha")
        theta_histogram(c_theta, f"{graphs}/c_theta")

        # relationship between alpha and theta got from all observations
        correlation_plot(np.concatenate([h_alpha, e_alpha, c_alpha]), np.concatenate([h_theta, e_theta, c_theta]), f"{graphs}/correlation")
