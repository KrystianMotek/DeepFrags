import argparse
import warnings
import os
import random

warnings.filterwarnings("ignore")
os.environ["TF_CPP_MIN_LOG_LEVEL"] = "3"

import tensorflow as tf
import numpy as np
from data.data import Label
from statistics.utils import Output

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-aa", type=str, help="amino acids sequence")
    parser.add_argument("-ss", type=str, help="secondary structure")
    parser.add_argument("-m", "--model", type=str, help="model to be used")
    args = parser.parse_args()

    np.set_printoptions(formatter={"float": lambda x: "{0:0.2f}".format(x)})
    tf.compat.v1.logging.set_verbosity(tf.compat.v1.logging.ERROR)

    aa = args.aa 
    ss = args.ss 
    model = args.model 

    label_vector = Label(aa, ss).to_vector() # generate vector for given sequence and secondary structure

    label_dim = len(label_vector) # number of features in label vector

    label = tf.keras.layers.Input(shape=(label_dim,)) # define input layer

    # latent space processing
    latent_variables = np.load(f"{model}/latent.npy")
    latent_dim = np.shape(latent_variables)[2]
    mean, log_variance = latent_variables[0], latent_variables[1]
    z = np.concatenate([i + np.exp(0.5 * j) * np.random.normal(loc=0.0, scale=1.0, size=(1, latent_dim)) for i, j in zip(mean, log_variance)])

    latent_sample = np.array(random.sample(list(z), 1)) # random sample from the latent space

    decoder = tf.keras.models.load_model(f"{model}/decoder.h5") # load decoder 

    decoder.load_weights(f"{model}/weights.h5", by_name=True, skip_mismatch=True) # take values of trained weights 

    output_vector = decoder.predict(tf.keras.layers.concatenate([latent_sample, label_vector]))[0] # predict cooridnates for given label

    output = Output(output_vector).to_original() # reverse output to original form
    
    displacement_norm = np.linalg.norm(output[0:3]) # compute norm of displacement

    # show results
    print(output)
    print("%.2f" % displacement_norm)
