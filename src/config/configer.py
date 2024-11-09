import os
import torch
import torch.nn as nn
from model import model_list
import logging

logging.basicConfig(level=logging.INFO)  # Set the logging level to INFO


def build_model(args):
    device = torch.device(
        "cuda" if torch.cuda.is_available() else "cpu"
    )
    model = model_list.RNN(args)

    network = model.to(device)
    total_params = network.get_num_params()

    logging.info(f"Network loaded to device {device}, device num {torch.cuda.device_count()}")
    logging.info(f"Total number of parameters: {total_params}")

    return network



