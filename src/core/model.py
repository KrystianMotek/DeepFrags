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
    return 0.5 * tf.reduce_mean(tf.square(mean) + tf.exp(log_variance) - log_variance - 1)


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
        self.label_dim = (20 + 3) * self.n + 3

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
        self.read_parameters()
        self.load_data()

        self.model = CVAE(n=self.n, latent_dim=self.latent_dim, encoder_h=self.encoder_h, decoder_h=self.decoder_h)

        logging.info(f"Trainer initiated from {self.config}")

    def read_parameters(self):
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
        self.beta = parameters["beta"]

        config_file.close()

    def load_data(self):
        inputs = np.load(self.inputs)
        labels = np.load(self.labels)
        
        self.training_inputs = inputs[:self.observations]
        self.training_labels = labels[:self.observations]

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
            indices = tf.random.shuffle(tf.range(0, self.observations, dtype=tf.int32))
            training_inputs = tf.gather(self.training_inputs, indices)
            training_labels = tf.gather(self.training_labels, indices)

            optimizer = tf.keras.optimizers.Adam(learning_rate=self.learning_rate)

            steps = int(self.observations / self.batch)
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

            training_message = f"total {total:.6f} reconstruction {reconstruction:.6f} kl {kl:.6f}"
            
            logging.info(f"Epoch {epoch+1}/{self.epochs}" + " " + training_message)
        
        logging.info("Fitting completed")

    def save(self):
        work_directory = os.path.dirname(self.config)
        
        self.model.encoder.save(f"{work_directory}/encoder.pb")
        self.model.decoder.save(f"{work_directory}/decoder.pb")
        
        self.model.encoder.save(f"{work_directory}/encoder.h5")
        self.model.decoder.save(f"{work_directory}/decoder.h5")
        
        self.model.encoder.save_weights(f"{work_directory}/encoder_weights.h5")
        self.model.decoder.save_weights(f"{work_directory}/decoder_weights.h5")

        # latent space variables
        np.save(f"{work_directory}/latent.npy", self.model.encode(self.training_inputs, self.training_labels))
    

class DecoderLoader:
    def __init__(self, decoder, latent):
        self.decoder = decoder
        self.latent = latent
        
        self.decoder = tf.keras.models.load_model(self.decoder)

        self.latent = np.load(self.latent)[0] # load samples from the latent space

    def predict(self, labels):
        z = np.array(random.choices(list(self.latent), k=len(labels)))
        return self.decoder.predict(tf.keras.layers.concatenate([z, labels]))