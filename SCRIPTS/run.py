import argparse
import warnings
import os
import random

warnings.filterwarnings("ignore")
os.environ["TF_CPP_MIN_LOG_LEVEL"] = "3"

import tensorflow as tf
import numpy as np
from statistics.graphs import *
from statistics.statistics import *

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-m", "--model", type=str, help="model to be used")
    parser.add_argument("-l", "--labels", type=str, help="sequences and secondary structures saved in a binary format")
    args = parser.parse_args()

    np.set_printoptions(formatter={"float": lambda x: "{0:0.2f}".format(x)})
    tf.compat.v1.logging.set_verbosity(tf.compat.v1.logging.ERROR)

    model = args.model 
    model_name = os.path.splitext(os.path.basename(model))[0]
    labels = np.load(args.labels)

    # labels dimensionality
    observations = np.shape(labels)[0] 
    label_dim = np.shape(labels)[1] 

    label = tf.keras.layers.Input(shape=(label_dim,)) # define input layer

    # latent space processing
    latent_variables = np.load(f"{model}/latent.npy")
    latent_dim = np.shape(latent_variables)[2]
    mean, log_variance = latent_variables[0], latent_variables[1]
    z = np.concatenate([i + np.exp(0.5 * j) * np.random.normal(loc=0.0, scale=1.0, size=(1, latent_dim)) for i, j in zip(mean, log_variance)])

    # sample from the latent space
    latent_sample = np.array(random.sample(list(z), observations)) 

    # load decoder with trained parameters
    decoder = tf.keras.models.load_model(f"{model}/decoder.h5")
    decoder.load_weights(f"{model}/weights.h5", by_name=True, skip_mismatch=True)

    # reconstruct cooridnates 
    raw_output = decoder.predict(tf.keras.layers.concatenate([latent_sample, labels]))
    outputs = [Output(vector) for vector in raw_output]

    # get secondary structures from labels 
    ss = [extract_ss(vector) for vector in labels] 

    # secondary structure connected with plane and dihedral angles
    angles = np.concatenate([angles_distribution(ss, output) for ss, output in zip(ss, outputs)]) 

    h_angles, e_angles, c_angles = hec_distribution(angles)[0], hec_distribution(angles)[1], hec_distribution(angles)[2]

    h_alpha = [float(angle[0]) for angle in h_angles]
    h_theta = [float(angle[1]) for angle in h_angles]

    e_alpha = [float(angle[0]) for angle in e_angles]
    e_theta = [float(angle[1]) for angle in e_angles]

    c_alpha = [float(angle[0]) for angle in c_angles]
    c_theta = [float(angle[1]) for angle in c_angles]

    # generate plots and save them in created directory

    os.makedirs(f"{model}/graphs")

    alpha_histogram(h_alpha, f"{model}/graphs//h_alpha")
    theta_histogram(h_theta, f"{model}/graphs/h_theta")
    
    alpha_histogram(e_alpha, f"{model}/graphs/e_alpha")
    theta_histogram(e_theta, f"{model}/graphs/e_theta")

    alpha_histogram(c_alpha, f"{model}/graphs/c_alpha")
    theta_histogram(c_theta, f"{model}/graphs/c_theta")

    correlation_plot(np.concatenate([h_alpha, e_alpha, c_alpha]), np.concatenate([h_theta, e_theta, c_theta]), f"{model}/graphs/correlation")
