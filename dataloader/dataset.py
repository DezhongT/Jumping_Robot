import random

import torch
from numpy.random import normal as gen_normal
from os import path as osp
import numpy as np
import logging
from torch.utils.data import Dataset, DataLoader
import math

class RobotData(Dataset):
    def __init__(self, np_data, mean, std):
        self.data = []
        for i in range(np_data.shape[0]):
            alpha = (np_data[i, 0] - mean[0])/std[0]
            # alpha = np_data[i, 0]
            # alpha *= (180/ math.pi)
            # alpha = alpha.astype(np.float32)
            # alpha = np_data[i, 0]
            compressL = (np_data[i, 1] - mean[1])/std[1]
            # compressL = np_data[i, 1]
            mu = (np_data[i, 2] - mean[2])/std[2]
            H = (np_data[i, 3] - mean[3])/std[3]
            rho = (np_data[i, 4] - mean[4])/std[4]
            L2 = (np_data[i, 5] - mean[5])/std[5]
            # L2 = np_data[i, 5]
            y = (np_data[i, 6] - mean[6])/std[6]
            x = (np_data[i, 7] - mean[7])/std[7]

            print(i, alpha)
            # features = [x, y, mu, rho]

            features = [x, y]
            output = [alpha, compressL]
            # output = [alpha, compressL, L2]
            self.data.append((torch.Tensor(features), torch.tensor(output)))

    def __len__(self):
        return len(self.data)

    def __getitem__(self, idx):
        return self.data[idx]


