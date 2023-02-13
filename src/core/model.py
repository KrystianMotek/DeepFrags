import os
import json
import random
import logging
import warnings

warnings.filterwarnings("ignore")
os.environ['TF_CPP_MIN_LOG_LEVEL'] = '3' 

import tensorflow as tf
import numpy as np


def kl_loss(mean, log_variance):
    return tf.reduce_mean(0.5 * (tf.square(mean) + tf.exp(log_variance) - log_variance - 1))


def reconstruction_loss(x, x_reconstructed):
    return tf.reduce_mean(tf.keras.losses.mean_squared_error(x, x_reconstructed))


def latent_sample(mean, log_variance):
    epsilon = tf.random.normal(shape=(1,))
    return mean + tf.exp(0.5 * log_variance) * epsilon


class CVAE(tf.keras.Model):
    def __init__(self, n, encoder_h, decoder_h, latent_dim):
        super().__init__()
        self.n = n
        self.encoder_h = encoder_h
        self.decoder_h = decoder_h
        self.latent_dim = latent_dim

        # expected data dimensions
        self.input_dim = 3 * self.n
        self.label_dim = (20 + 3) * self.n + 1

        self.encoder = self.build_encoder()
        self.decoder = self.build_decoder()

    def encode(self, inputs, labels):
        data = tf.keras.layers.concatenate([inputs, labels])
        mean, log_variance = tf.split(self.encoder(data), num_or_size_splits=2, axis=1)
        return mean, log_variance

    def decode(self, z, labels):
        data = tf.keras.layers.concatenate([z, labels])
        return self.decoder(data)

    def build_encoder(self):
        encoder = tf.keras.Sequential(name="encoder")
        encoder.add(tf.keras.layers.InputLayer([self.input_dim + self.label_dim]))
        encoder.add(tf.keras.layers.Dense(self.encoder_h, activation="relu"))
        encoder.add(tf.keras.layers.Dense(2 * self.latent_dim, activation="linear"))
        return encoder

    def build_decoder(self):
        decoder = tf.keras.Sequential(name="decoder")
        decoder.add(tf.keras.layers.InputLayer([self.latent_dim + self.label_dim]))
        decoder.add(tf.keras.layers.Dense(self.decoder_h, activation="relu"))
        decoder.add(tf.keras.layers.Dense(self.input_dim, activation="linear"))
        return decoder


class Trainer:
    def __init__(self, config):
        self.config = config
        self.load_parameters()
        self.split_data()

        # create model
        self.model = CVAE(n=self.n, latent_dim=self.latent_dim, encoder_h=self.encoder_h, decoder_h=self.decoder_h)

        logging.info(f"Trainer initiated from {self.config}")

    def load_parameters(self):
        config_file = open(self.config, "r")
        parameters = json.loads(config_file.read())
    
        self.n = parameters["n"]
        self.encoder_h = parameters["encoder_h"]
        self.decoder_h = parameters["decoder_h"]
        self.latent_dim = parameters["latent_dim"]
        self.observations = parameters["observations"]
        self.learning_rate = parameters["learning_rate"]
        self.inputs = parameters["inputs"]
        self.labels = parameters["labels"]
        self.epochs = parameters["epochs"]
        self.batch = parameters["batch"]
        self.split = parameters["split"]
        self.beta = parameters["beta"]

        config_file.close()

    def split_data(self):
        # divide data into training and validation sets
        inputs = np.load(self.inputs)
        labels = np.load(self.labels)

        self.split_index = int(self.split * self.observations)
        
        self.training_inputs = inputs[0:self.split_index]
        self.training_labels = labels[0:self.split_index]

        self.validation_inputs = inputs[self.split_index:]
        self.validation_labels = labels[self.split_index:]

    def losses(self, inputs, labels):
        # pass data through the network
        mean, log_variance = self.model.encode(inputs, labels)
        z = latent_sample(mean, log_variance)
        inputs_reconstructed = self.model.decode(z, labels)

        # compute losses
        kl = kl_loss(mean, log_variance)
        reconstruction = reconstruction_loss(inputs, inputs_reconstructed)
        total = self.beta * kl + reconstruction

        return total, reconstruction, kl

    def train(self):
        for epoch in range(self.epochs):
            # data shuffling before each epoch
            indices = tf.random.shuffle(tf.range(0, self.split_index, dtype=tf.int32))
            training_inputs = tf.gather(self.training_inputs, indices)
            training_labels = tf.gather(self.training_labels, indices)

            optimizer = tf.keras.optimizers.Adam(learning_rate=self.learning_rate)

            steps = int(self.split_index / self.batch)
            for step in range(steps):
                start = int(step * self.batch)
                end = int(start + self.batch)
                batch_inputs = training_inputs[start:end]
                batch_labels = training_labels[start:end]

                with tf.GradientTape() as tape:
                    total, reconstruction, kl = self.losses(batch_inputs, batch_labels)

                # update trainable parameters
                grads = tape.gradient(total, self.model.trainable_weights)
                optimizer.apply_gradients(zip(grads, self.model.trainable_weights))

            val_total, val_reconstruction, val_kl = self.losses(self.validation_inputs, self.validation_labels)

            training_message = f"total {total:.4f} reconstruction {reconstruction:.4f} kl {kl:.4f}"
            validation_message = f"val_total {val_total:.4f} val_reconstruction {val_reconstruction:.4f} val_kl {val_kl:.4f}"
            
            logging.info(f"Epoch {epoch+1}/{self.epochs}" + " " + training_message + " " + validation_message)
        
        logging.info("Fitting completed")

    def save(self):
        work_directory = os.path.dirname(self.config)

        # models and obtained weights
        self.model.encoder.save(f"{work_directory}/encoder.h5")
        self.model.decoder.save(f"{work_directory}/decoder.h5")
        self.model.save_weights(f"{work_directory}/weights.h5")

        # latent space variables
        np.save(f"{work_directory}/latent.npy", self.model.encode(self.training_inputs, self.training_labels))
    

class DecoderLoader:
    def __init__(self, decoder, weights, latent):
        self.decoder = decoder
        self.weights = weights
        self.latent = latent
        
        self.decoder = tf.keras.models.load_model(self.decoder)
        self.decoder.load_weights(self.weights, by_name=True, skip_mismatch=True)

        self.latent_processing() # load samples from the latent space

    def latent_processing(self):
        latent_variables = np.load(self.latent)
        mean, log_variance = latent_variables[0], latent_variables[1]
        self.z = latent_sample(mean, log_variance)

    def predict(self, labels):
        z = np.array(random.sample(list(self.z), len(labels)))
        return self.decoder.predict(tf.keras.layers.concatenate([z, labels]))