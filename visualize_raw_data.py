import numpy as np
import matplotlib.pyplot as plt
import yaml
# from mpl_toolkits.mplot3d import Axes3D
from config import configer
import os
import logging
from dataloader import dataset
from sklearn.model_selection import train_test_split
from torch.utils.data import Dataset, DataLoader
import torch.nn as nn
import torch
import torch.optim as optim
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

from os import path as osp


np_data = np.loadtxt("data.txt")
mean = np.mean(np_data, axis = 0)
std = np.std(np_data, axis = 0)

alpha = (np_data[:, 0] - mean[0])/std[0]
compressL = (np_data[:, 1] - mean[1])/std[1]
mu = (np_data[:, 2] - mean[2])/std[2]
H = (np_data[:, 3] - mean[3])/std[3]
rho = (np_data[:, 4] - mean[4])/std[4]
L2 = (np_data[:, 5] - mean[5])/std[5]
y = (np_data[:, 6] - mean[6])/std[6]
x = (np_data[:, 7] - mean[7])/std[7]




fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
# ax.scatter(x[idx], y[idx], L2[idx], label="train")
ax.scatter(x, y, alpha, label="train")


ax.set_xlabel('X Label')
ax.set_ylabel('Y Label')
# ax.set_aspect('equal')

plt.legend()


plt.show()
