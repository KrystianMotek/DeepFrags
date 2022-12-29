import argparse
import warnings
import os
import random

warnings.filterwarnings("ignore")
os.environ["TF_CPP_MIN_LOG_LEVEL"] = "3"

import tensorflow as tf
import numpy as np
from data.data import Label
from statistics.statistics import Output

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-aa", type=str, help="amino acids sequence")
    parser.add_argument("-ss", type=str, help="secondary structure")
    parser.add_argument("-m", "--model", type=str, help="model to be used")
    args = parser.parse_args()

    tf.compat.v1.logging.set_verbosity(tf.compat.v1.logging.ERROR)

    aa = args.aa 
    ss = args.ss 
    model = args.model 

    label_vector = Label(aa, ss).to_vector() 

    label_dim = len(label_vector) # number of features in label vector

    label = tf.keras.layers.Input(shape=(label_dim,)) # define input layer

    # latent space processing
    latent_variables = np.load(f"{model}/latent.npy")
    latent_dim = np.shape(latent_variables)[2]
    mean, log_variance = latent_variables[0], latent_variables[1]
    z = np.concatenate([i + np.exp(0.5 * j) * np.random.normal(loc=0.0, scale=1.0, size=(1, latent_dim)) for i, j in zip(mean, log_variance)])

    # random sample from the latent space
    latent_sample = np.array(random.sample(list(z), 1)) 

    # load decoder 
    decoder = tf.keras.models.load_model(f"{model}/decoder.h5")  
    decoder.load_weights(f"{model}/weights.h5", by_name=True, skip_mismatch=True) 

    # predict cooridnates for given label
    output_vector = decoder.predict(tf.keras.layers.concatenate([latent_sample, label_vector]))[0] 

    # reverse output to PDB format
    output = Output(output_vector).to_pdb() 

    # show results
    for line in output:
        print(line)
