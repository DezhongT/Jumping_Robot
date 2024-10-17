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
from sklearn.neighbors import KernelDensity


def loadModel(yaml_path = "model.yaml", checkpoint_path = "./output/checkpoints/inverse_model.pt"):
    with open(yaml_path, 'r') as file:
        args = yaml.safe_load(file)
    model = configer.build_model(args)
    checkpoint = torch.load(checkpoint_path, weights_only=True)
    model.load_state_dict(checkpoint['model_state_dict'])

    return model


global mean_alpha, std_alpha, mean_compressL, std_compressL, \
mean_mu, std_mu, mean_rho, std_rho, mean_x, std_x, mean_y, std_y \

def loadData(data_path = "train_data.txt"):
    np_data = np.loadtxt(data_path)
    mean = np.mean(np_data, axis = 0)
    std = np.std(np_data, axis = 0)

    alpha = np_data[:, 0]
    compressL = np_data[:, 1]
    mu = np_data[:, 2]
    # H = np_data[:, 3]
    rho = np_data[:, 4]
    # L2 = np_data[:, 5]
    y = np_data[:, 6]
    x  = np_data[:, 7]

    global mean_alpha, std_alpha, mean_compressL, std_compressL, \
        mean_mu, std_mu, mean_rho, std_rho, mean_x, std_x, mean_y, std_y

    mean_alpha = mean[0]
    std_alpha = std[0]
    mean_compressL = mean[1]
    std_compressL = std[1]
    mean_mu = mean[2]
    std_mu = std[2]
    mean_rho = mean[4]
    std_rho = std[4]
    mean_y = mean[6]
    std_y = std[6]
    mean_x = mean[7]
    std_x = std[7]

    return x, y, mu, rho, alpha, compressL

def printInfo():
    global mean_alpha, std_alpha, mean_compressL, std_compressL, \
        mean_mu, std_mu, mean_rho, std_rho, mean_x, std_x, mean_y, std_y

    print(mean_alpha, std_alpha, mean_compressL, std_compressL, \
          mean_mu, std_mu, mean_rho, std_rho, mean_y, std_y,  mean_x, mean_y)

def getRange(X):
    return np.amin(X), np.amax(X)


def generate_samples_within_bounds(kde, n_samples, lower_bound, upper_bound):
    samples = []
    while len(samples) < n_samples:
        new_samples = kde.sample(n_samples)
        within_bounds = np.all((new_samples >= lower_bound) & (new_samples <= upper_bound), axis=1)
        samples.append(new_samples[within_bounds])
    return np.vstack(samples)[:n_samples]


def visualize_samples(samples, input):
    fig = plt.figure()
    ax = fig.add_subplot(1, 2, 1, projection='3d')
    ax.scatter(input[:, 0], input[:, 1], input[:, 2], label="train data")
    ax.scatter(samples[:, 0], samples[:, 1], samples[:, 2], label="samples  ")
    ax.legend()

    ax1 = fig.add_subplot(1, 2, 2, projection='3d')
    ax1.scatter(input[:, 0], input[:, 1], input[:, 3], label="train data")
    ax1.scatter(samples[:, 0], samples[:, 1], samples[:, 3], label="samples  ")
    ax1.legend()
    plt.show()


def pred(model, features):
    model.eval()
    device = torch.device("cuda" if torch.cuda.is_available() else "cpu")

    global mean_alpha, std_alpha, mean_compressL, std_compressL, \
        mean_mu, std_mu, mean_rho, std_rho, mean_x, std_x, mean_y, std_y
    # alpoha, compressL, mu, rho
    features[:, 0] = (features[:, 0] - mean_alpha)/std_alpha
    features[:, 1] = (features[:, 1] - mean_compressL)/std_compressL
    features[:, 2] = (features[:, 2] - mean_mu)/std_mu
    features[:, 3] = (features[:, 3] - mean_rho)/std_rho

    features = features.astype(np.float32)
    features = torch.tensor(features).to(device)

    output = model(features)
    output = output.detach().cpu().numpy()
    output *= 0.02

    return output


from scipy.spatial import Delaunay

def checkTarget(points, target, plot_flag = False):
    tri = Delaunay(points)
    edges = set()
    def add_edge(i, j):
        """ Add a sorted edge, i < j """
        if (i, j) in edges or (j, i) in edges:
            if (j, i) in edges:
                edges.remove((j, i))
            return
        edges.add((i, j))

    # Loop over the triangles
    for ia, ib, ic in tri.simplices:
        add_edge(ia, ib)
        add_edge(ib, ic)
        add_edge(ic, ia)
    from collections import defaultdict
    def extract_contours(edges):
        neighbors = defaultdict(set)
        for i, j in edges:
            neighbors[i].add(j)
            neighbors[j].add(i)

        def find_contour(start):
            contour = [start]
            current = start
            previous = None
            while True:
                next_point = (neighbors[current] - {previous}).pop()
                if next_point == start:
                    break
                contour.append(next_point)
                previous, current = current, next_point
            return contour

        contours = []
        used_points = set()
        for point in neighbors:
            if point not in used_points:
                contour = find_contour(point)
                contours.append(contour)
                used_points.update(contour)
        return contours

    contours = extract_contours(edges)

    contour = contours[0]
    vertices = points[contour, :2]
    flag = is_point_in_polygon(target.reshape(-1, ), vertices)

    if plot_flag:
        for contour in contours:
            plt.fill(points[contour, 0], points[contour, 1], alpha=0.5, label = "Available region")

    return flag



def is_point_in_polygon(point, vertices):
    x, y = point
    n = len(vertices)
    inside = False

    p1x, p1y = vertices[0]
    for i in range(n + 1):
        p2x, p2y = vertices[i % n]
        if y > min(p1y, p2y):
            if y <= max(p1y, p2y):
                if x <= max(p1x, p2x):
                    if p1y != p2y:
                        xinters = (y - p1y) * (p2x - p1x) / (p2y - p1y) + p1x
                    if p1x == p2x or x <= xinters:
                        inside = not inside
        p1x, p1y = p2x, p2y

    return inside


def inverseOptimize(inputs, model, target_p, max_iter = 1000, verbose = False):
    global mean_alpha, std_alpha, mean_compressL, std_compressL, \
        mean_mu, std_mu, mean_rho, std_rho, mean_x, std_x, mean_y, std_y
    device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
    loss_fn = torch.nn.MSELoss()
    target_p = target_p.astype(np.float32)
    target_p = target_p/0.02
    target = torch.tensor(target_p).view(-1, 2 ).to(device)

    inputs[0] = (inputs[0] - mean_alpha)/std_alpha
    inputs[1] = (inputs[1] - mean_compressL)/std_compressL
    inputs[2] = (inputs[2] - mean_mu)/std_mu
    inputs[3] = (inputs[3] - mean_rho)/std_rho

    constant_inputs = torch.tensor([inputs[2], inputs[3]])
    variable_inputs = torch.tensor([inputs[0], inputs[1]], requires_grad=True)
    optimizer = optim.Adam([variable_inputs], lr=0.001)  # lr is the learning rate1


    for i in range(max_iter):
        inputs= torch.cat((variable_inputs, constant_inputs)).unsqueeze(0).to(device)
        inputs = inputs.to(torch.float32)
        outputs = model(inputs)
        loss = loss_fn(target, outputs)

        optimizer.zero_grad()
        loss.backward()
        optimizer.step()
        if verbose and (i + 1) % 100 == 0:
            print(f"Iter {i+1}: Loss is {loss}")

    print(f"Inverse design optimization loss is {loss}")
    constant_inputs = constant_inputs.detach().cpu().numpy()
    variable_inputs = variable_inputs.detach().cpu().numpy()

    mu = constant_inputs[0] * std_mu + mean_mu
    rho = constant_inputs[1] * std_rho + mean_rho
    alpha = variable_inputs[0] * std_alpha + mean_alpha
    compressL = variable_inputs[1] * std_compressL + mean_compressL

    return alpha, compressL, mu, rho

def data_driven(test_samples):
    x, y, mu, rho, alpha, compressL = loadData()
    alpha_f, alpha_c = getRange(alpha)
    compressL_f, compressL_c = getRange(compressL)
    rho_f, rho_c = getRange(rho)
    mu_f, mu_c = getRange(mu)

    # load the data-driven model
    yaml_path = "model.yaml"
    checkpoint_path = "./output/checkpoints/forward_model.pt"
    model = loadModel(yaml_path=yaml_path, checkpoint_path=checkpoint_path)

    for sample in test_samples:
        print(sample)
        x, y, mu, rho = sample
    exit(0)
