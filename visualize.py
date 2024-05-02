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


yaml_path = "model.yaml"
with open(yaml_path, 'r') as file:
    args = yaml.safe_load(file)

model = configer.build_model(args)
checkpoint_path = "./output/checkpoints/checkpoint_200.pt"
checkpoint = torch.load(checkpoint_path)
model.load_state_dict(checkpoint['model_state_dict'])

np_data = np.loadtxt("data.txt")
mean = np.mean(np_data, axis = 0)
std = np.std(np_data, axis = 0)




pred = np.zeros_like(np_data)
# mean = np.zeros_like(mean)
# std = np.ones_like(std)


model.eval()
alpha = (np_data[:, 0] - mean[0])/std[0]
compressL = (np_data[:, 1] - mean[1])/std[1]
mu = (np_data[:, 2] - mean[2])/std[2]
H = (np_data[:, 3] - mean[3])/std[3]
rho = (np_data[:, 4] - mean[4])/std[4]
L2 = (np_data[:, 5] - mean[5])/std[5]
y = (np_data[:, 6] - mean[6])/std[6]
x = (np_data[:, 7] - mean[7])/std[7]

# features = np.hstack([x.reshape(-1, 1), y.reshape(-1, 1), mu.reshape(-1, 1)])
features = np.hstack([x.reshape(-1, 1), y.reshape(-1, 1)])
features = features.astype(np.float32)
device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
features = torch.tensor(features).to(device)
pred = model(features)
pred = pred.detach().cpu().numpy()


pred[:, 0] = pred[:, 0] * std[0] + mean[0]
pred[:, 1] = pred[:, 1] * std[1] + mean[1]
output = np.zeros_like(np_data)
output[:, 0] = np_data[:, 0]
output[:, 1] = np_data[:, 1]




mu = np_data[:, 2]
unique_mu = np.unique(mu)
unique_rho = np.unique(rho)
unique_compressL = np.unique(compressL)
unique_L2 = np.unique(L2)
unique_H = np.unique(H)
unique_alpha = np.unique(alpha)

print(unique_mu)
idx = (mu == unique_mu[0])

fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
ax.scatter(x, y, output[:, 0], label="train")
ax.scatter(x, y, pred[:, 0], label="pred")

output = np.hstack([alpha.reshape(-1, 2)])

ax.set_xlabel('X Label')
ax.set_ylabel('Y Label')
ax.set_zlabel("alpha")
# ax.set_aspect('equal')
ax.view_init(elev=30, azim=50)  # Change the elevation and azimuth angles here

plt.legend()
plt.savefig('alpha.png')  # You can specify the file format here, like .png, .jpg, .pdf, etc.


plt.show()
