import os
import argparse
import logging
from core.model import Trainer

logging.getLogger("tensorflow").disabled=True
logging.getLogger("h5py._conv").disabled=True

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-cfg", "--config", type=str, help="configuration file path")
    args = parser.parse_args()

    logging_file = f"{os.path.dirname(args.config)}/logging"
    logging_format = "%(asctime)s %(name)s %(message)s"

    logging.basicConfig(filename=logging_file, filemode="a", format=logging_format, datefmt="%Y-%m-%d %H:%M:%S", level=logging.DEBUG)

    trainer = Trainer(config=args.config)

    trainer.train()
    trainer.save()