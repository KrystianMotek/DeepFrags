import argparse
import logging
import numpy as np
from model import DecoderLoader
from features import LabelMLP
from statistical import Output

logging.getLogger("tensorflow").disabled=True
logging.getLogger("h5py._conv").disabled=True

if __name__ == "__main__":
    pass