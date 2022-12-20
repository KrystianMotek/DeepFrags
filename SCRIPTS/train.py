import argparse
import json
import warnings
import os
import logging

warnings.filterwarnings("ignore")
os.environ["TF_CPP_MIN_LOG_LEVEL"] = "3"

import tensorflow as tf
import numpy as np


def latent_sample(mean, log_variance):
    epsilon = tf.random.normal(shape=(1,))
    return mean + tf.exp(0.5 * log_variance) * epsilon


def kl_loss(mean, log_variance):
    return tf.reduce_mean(0.5 * (tf.square(mean) + tf.exp(log_variance) - log_variance - 1))


def mse_loss(x, x_reconstructed):
    return tf.reduce_mean(tf.keras.losses.mean_squared_error(x, x_reconstructed))


if __name__ == "__main__":
    tf.compat.v1.disable_eager_execution()

    parser = argparse.ArgumentParser()
    parser.add_argument("-cfg", "--config", type=str, help="configuration file path")
    args = parser.parse_args()

    # read configuration file and load parameters
    config = args.config
    config_file = open(config, "r")
    parameters = json.loads(config_file.read())
    config_file.close()

    batch = parameters["batch"]
    epochs = parameters["epochs"]
    observations = parameters["observations"]
    learning_rate = parameters["learning_rate"]

    # number of steps in a single epoch
    steps = int(observations / batch)

    # latent space dimension
    latent_dim = parameters["latent_dim"]

    # scales for each term in the loss function
    kl_factor = parameters["kl_factor"]
    mse_factor = parameters["mse_factor"]

    # nodes in hidden layers
    encoder_hidden = parameters["encoder_hidden"]
    decoder_hidden = parameters["decoder_hidden"]

    # load training data
    inputs = np.load(parameters["inputs"])
    labels = np.load(parameters["labels"])

    # number of features in each tensor
    x_dim = np.shape(inputs)[1]
    label_dim = np.shape(labels)[1]

    # directory dedicated to work files
    work_directory = os.path.dirname(config)
    
    # logging configuration
    logging_file = f"{work_directory}/logging"
    logging_format = "%(asctime)s %(name)s %(message)s"
    logging.getLogger("tensorflow").disabled=True
    logging.getLogger("h5py._conv").disabled=True
    logging.basicConfig(filename=logging_file, filemode="a", format=logging_format, datefmt="%Y-%m-%d %H:%M:%S", level=logging.DEBUG)

    # define input layers
    x = tf.keras.layers.Input(shape=(x_dim,))
    label = tf.keras.layers.Input(shape=(label_dim,))

    # encoder architecture
    encoder = tf.keras.Sequential(name="encoder")
    encoder.add(tf.keras.layers.InputLayer([x_dim + label_dim]))
    encoder.add(tf.keras.layers.Dense(encoder_hidden, activation="sigmoid"))
    encoder.add(tf.keras.layers.Dense(2 * latent_dim, activation="linear"))

    # decoder architecture
    decoder = tf.keras.Sequential(name="decoder")
    decoder.add(tf.keras.layers.InputLayer([latent_dim + label_dim]))
    decoder.add(tf.keras.layers.Dense(decoder_hidden, activation="sigmoid"))
    decoder.add(tf.keras.layers.Dense(x_dim, activation="linear"))

    # split encoded data into two vectors
    mean, log_variance = tf.split(encoder(tf.keras.layers.concatenate([x, label])), num_or_size_splits=2, axis=1)

    # sample from the latent space
    z = latent_sample(mean, log_variance)

    # predict original data
    x_reconstructed = decoder(tf.keras.layers.concatenate([z, label]))

    def total_loss(mean, log_variance):
        kl = kl_loss(mean, log_variance) 
        mse = mse_loss(x, x_reconstructed)
        return kl_factor * kl + mse_factor * mse
    
    # build model
    cvae = tf.keras.Model([x, label], x_reconstructed)

    logging.info(f"CVAE training initiated from {os.path.splitext(os.path.basename(config))[0]}.cfg")

    # compile and fit
    cvae.compile(optimizer=tf.keras.optimizers.Adam(learning_rate=learning_rate), loss=total_loss, metrics=[kl_loss, mse_loss])
    history = cvae.fit([inputs, labels], inputs, batch_size=batch, epochs=epochs, verbose=2, shuffle=True, steps_per_epoch=steps)

    logging.info(f"Training completed")

    # save metrics
    metrics = np.transpose([history.history["kl_loss"], history.history["mse_loss"]])
    np.save(f"{work_directory}/metrics", metrics)

    latent_variables = encoder.predict(tf.keras.layers.concatenate([inputs, labels]), steps=1)

    # save variables from the latent space
    np.save(f"{work_directory}/latent", np.split(latent_variables, indices_or_sections=2, axis=1))

    # save architectures and obtained weights
    encoder.save(f"{work_directory}/encoder.h5")
    decoder.save(f"{work_directory}/decoder.h5")
    cvae.save_weights(f"{work_directory}/weights.h5")
