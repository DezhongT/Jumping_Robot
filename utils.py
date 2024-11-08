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
import tracemalloc
import time
import subprocess

# bayesian opt package
from argparse import Namespace
from dragonfly import load_config
from dragonfly.exd.experiment_caller import CPFunctionCaller
from dragonfly.opt import gp_bandit


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
    if verbose:
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

    tracemalloc.start()
    start_time = time.time()

    # load the data-driven model
    yaml_path = "model.yaml"
    checkpoint_path = "./output/checkpoints/forward_model.pt"
    model = loadModel(yaml_path=yaml_path, checkpoint_path=checkpoint_path)

    cache = []
    for sample in test_samples:
        mu_t, rho_t, x_t, y_t = sample
        # check the available region
        alpha = np.linspace(alpha_f, alpha_c, 10)
        compressL = np.linspace(compressL_f, compressL_c, 10)

        alpha, compressL= np.meshgrid(alpha, compressL, indexing='ij')
        alpha_input = alpha.flatten()
        compressL_input = compressL.flatten()
        mu_input = np.ones_like(alpha_input) * mu_t
        rho_input = np.ones_like(alpha_input) * rho_t

        features = np.vstack([alpha_input, compressL_input, mu_input, rho_input]).T
        outputs = pred(model, features)
        target_p = np.asarray([x_t, y_t]).reshape(-1, 2)
        flag = checkTarget(outputs, target_p, plot_flag=False)
        if not flag:
            continue

        # select initial guess
        offset = outputs - target_p
        offset = np.sqrt((offset**2).sum(axis = 1))

        min_index = np.argmin(offset)
        alpha_int = alpha_input[min_index]
        compressL_int = compressL_input[min_index]

        inputs = [alpha_int, compressL_int, mu_t, rho_t]
        alpha_r, compressL_r, mu_r, rho_r = inverseOptimize(inputs, model, target_p.copy())
        cache.append([alpha_r, compressL_r, mu_r, rho_r, target_p])
    current, peak = tracemalloc.get_traced_memory()
    tracemalloc.stop()
    peak_memory = peak/(1024 ** 2)
    elapsed_time = time.time() - start_time

    # evaluate the accuracy
    accuracy = 0
    for data in cache:
        alpha_r, compressL_r, mu_r, rho_r, target_p = data
        cmd = './simulations/simDER ./simulations/option.txt'
        suffix = f" -- mu {mu_r} -- totalMass {rho_r} -- angleRight {alpha_r} -- compressRatio {compressL_r}"
        cmd = cmd + suffix
        BASimProcess = subprocess.Popen(cmd, shell=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
        BASimProcess.communicate()
        fileName = f"datafiles/simDiscreteNet_l2_0.010000_compressRatio_{compressL_r:.6f}_angleRight_{alpha_r:.6f}_mu_{mu_r:.6f}_mass_{rho_r:.6f}.txt"
        data1 = np.loadtxt(fileName)
        idx = np.argmax(data1[:, 2])
        node= data1[idx]
        y = np.linalg.norm(target_p - node[1:])
        accuracy +=  y
    results = {"accuracy": accuracy/len(test_samples), "memory_usage": peak_memory, "time_cost": elapsed_time/len(test_samples)}

    return results

def bayesian_opt(test_samples):
    x, y, mu, rho, alpha, compressL = loadData()
    alpha_f, alpha_c = getRange(alpha)
    compressL_f, compressL_c = getRange(compressL)

    tracemalloc.start()
    start_time = time.time()
    num_init = 6
    batch_size = 10

    options = Namespace(
        build_new_model_every = batch_size,
        init_captial = num_init-1,
        gpb_hp_tune_criterion = "ml-post_sampling",
    )

    domain_vars = [{'type': 'float',  'min': alpha_f, 'max': alpha_c, 'dim': 1}, # alpha
                   {'type': 'float',  'min': compressL_f, 'max': compressL_c, 'dim': 1}, # compressL
                   ]
    config_params = {"domain": domain_vars}
    config = load_config(config_params)

    print("Executing bayesian optimization...")
    cache_result = []
    accuracy = 0
    for sample in test_samples:
        mu_t, rho_t, x_t, y_t = sample
        target_p = np.asarray([x_t, y_t]).reshape(-1, 2)

        cache = {}

        func_caller = CPFunctionCaller(None, config.domain, domain_orderings=config.domain_orderings)
        opt = gp_bandit.CPGPBandit(func_caller, 'default', ask_tell_mode=True, options=options)  # opt is the optimizer object
        opt.initialise()
        initial_samples = opt.ask(num_init)
        for sample in initial_samples:
            alpha_t = sample[0][0]
            compressL_t = sample[1][0]
            # check if the file exist
            if (alpha_t, compressL_t) in cache:
                opt.tell([(sample, -cache[(alpha_t, compressL_t)])])
            else:
                cmd = './simulations/simDER ./simulations/option.txt'
                suffix = f" -- mu {mu_t} -- totalMass {rho_t} -- angleRight {alpha_t} -- compressRatio {compressL_t}"
                cmd = cmd + suffix
                BASimProcess = subprocess.Popen(cmd, shell=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
                BASimProcess.communicate()
                print(cmd)
                fileName = f"datafiles/simDiscreteNet_l2_0.010000_compressRatio_{compressL_t:.6f}_angleRight_{alpha_t:.6f}_mu_{mu_t:.6f}_mass_{rho_t:.6f}.txt"
                if os.path.getsize(fileName) > 0:
                    data1 = np.loadtxt(fileName)
                    idx = np.argmax(data1[:, 2])
                    node= data1[idx]
                    y = np.linalg.norm(target_p - node[1:])
                    opt.tell([(sample, -y)])
                    cache[(alpha_t, compressL_t)] = y
                os.remove(fileName)
        # update model
        opt._build_new_model()
        opt._set_next_gp()
        print("finished warm up")

        epochs = 0
        min_err = 1
        alpha_r = None
        compressL_r = None
        predicted = None
        while epochs < 10:
            # sampling
            batch_samples = []
            for i in range(batch_size):
                batch_samples.append(opt.ask())

            average_err = 0
            sample_num = 0
            for sample in batch_samples:
                alpha_t = sample[0][0]
                compressL_t = sample[1][0]
                # check if the file exist
                if (alpha_t, compressL_t) not in cache:
                    cmd = './simulations/simDER ./simulations/option.txt'
                    suffix = f" -- mu {mu_t} -- totalMass {rho_t} -- angleRight {alpha_t} -- compressRatio {compressL_t}"
                    cmd = cmd + suffix
                    BASimProcess = subprocess.Popen(cmd, shell=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
                    BASimProcess.communicate()
                    fileName = f"datafiles/simDiscreteNet_l2_0.010000_compressRatio_{compressL_t:.6f}_angleRight_{alpha_t:.6f}_mu_{mu_t:.6f}_mass_{rho_t:.6f}.txt"
                    if os.path.getsize(fileName) > 0:
                        data1 = np.loadtxt(fileName)
                        idx = np.argmax(data1[:, 2])
                        node= data1[idx]
                        y = np.linalg.norm(target_p - node[1:])
                        opt.tell([(sample, -y)])
                        cache[(alpha_t, compressL_t)] = y
                        average_err += cache[(alpha_t, compressL_t)]
                        sample_num += 1
                        if min_err > cache[(alpha_t, compressL_t)]:
                            min_err = cache[(alpha_t, compressL_t)]
                            alpha_r = alpha_t
                            compressL_r = compressL_t
                            predicted = node.copy()

                    os.remove(fileName)

            if sample_num > 0:
                average_err = average_err/sample_num
                if min_err < 3e-3 or average_err < 6e-3:
                    break
                print(f"epoch {epochs}, error: {average_err}, min_error : {min_err}, ({alpha_r} {compressL_r})")
                opt._build_new_model()
                opt._set_next_gp()
            epochs += 1

        print("finished evalute")
        cache_result.append([alpha_r, compressL_r, target_p, predicted])
        accuracy += min_err

    current, peak = tracemalloc.get_traced_memory()
    tracemalloc.stop()
    peak_memory = peak/(1024 ** 2)
    elapsed_time = time.time() - start_time
    
    results = {"accuracy": accuracy/len(test_samples), "memory_usage": peak_memory, "time_cost": elapsed_time/len(test_samples)}
    print(results)

    print(cache_result)
    return results


def gradient_descent(test_samples):
    tracemalloc.start()
    start_time = time.time()

    x, y, mu, rho, alpha, compressL = loadData()
    points = np.hstack((x.reshape(-1, 1), y.reshape(-1, 1), 100 * rho.reshape(-1, 1), mu.reshape(-1, 1)))
    accuracy = 0
    cache_result = []
    for sample in test_samples:
        mu_t, rho_t, x_t, y_t = sample

        target_p = np.asarray([x_t, y_t]).reshape(-1, 2)

        diff_p = points - np.asarray([x_t, y_t, 100 * rho_t, mu_t]).reshape(-1, 4)
        diff_p = np.sum(diff_p**2, 1)
        idx = np.argmin(diff_p)
        # initial guess
        alpha_g = alpha[idx]
        compressL_g = compressL[idx]
        print(idx, mu.shape, rho.shape, x.shape, y.shape)
        print(alpha_g, compressL_g, mu_t, rho_t, mu[idx], rho[idx], x[idx], y[idx], x_t, y_t)
        # exit(0)

        cmd = './simulations/simDER ./simulations/option.txt'
        suffix = f" -- mu {mu_t} -- totalMass {rho_t} -- angleRight {alpha_g} -- compressRatio {compressL_g}"
        cmd = cmd + suffix
        BASimProcess = subprocess.Popen(cmd, shell=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
        BASimProcess.communicate()
        fileName = f"datafiles/simDiscreteNet_l2_0.010000_compressRatio_{compressL_g:.6f}_angleRight_{alpha_g:.6f}_mu_{mu_t:.6f}_mass_{rho_t:.6f}.txt"
        data1 = np.loadtxt(fileName)
        idx = np.argmax(data1[:, 2])
        node= data1[idx]
        y_g = np.linalg.norm(target_p - node[1:])
        os.remove(fileName)
        print(y_g)

        iter = 0
        perturb = 1e-4

        while y_g > 1e-3 and iter < 100:
            # compute the numerical jacobian
            jacobian = np.zeros((2, ))
            suffix = f" -- mu {mu_t} -- totalMass {rho_t} -- angleRight {alpha_g + perturb} -- compressRatio {compressL_g}"
            cmd = cmd + suffix
            BASimProcess = subprocess.Popen(cmd, shell=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
            BASimProcess.communicate()
            fileName = f"datafiles/simDiscreteNet_l2_0.010000_compressRatio_{compressL_g:.6f}_angleRight_{alpha_g + perturb:.6f}_mu_{mu_t:.6f}_mass_{rho_t:.6f}.txt"
            data1 = np.loadtxt(fileName)
            idx = np.argmax(data1[:, 2])
            node= data1[idx]
            y_e = np.linalg.norm(target_p - node[1:])
            os.remove(fileName)

            jacobian[0] = (y_e - y_g)/perturb

            suffix = f" -- mu {mu_t} -- totalMass {rho_t} -- angleRight {alpha_g} -- compressRatio {compressL_g + perturb}"
            cmd = cmd + suffix
            BASimProcess = subprocess.Popen(cmd, shell=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
            BASimProcess.communicate()
            fileName = f"datafiles/simDiscreteNet_l2_0.010000_compressRatio_{compressL_g + perturb:.6f}_angleRight_{alpha_g:.6f}_mu_{mu_t:.6f}_mass_{rho_t:.6f}.txt"
            data1 = np.loadtxt(fileName)
            idx = np.argmax(data1[:, 2])
            node= data1[idx]
            y_e = np.linalg.norm(target_p - node[1:])
            os.remove(fileName)

            jacobian[1] = (y_e - y_g)/perturb

            if np.linalg.norm(jacobian) > 1e-2:
                jacobian = jacobian/np.linalg.norm(jacobian) * 1e-2

            step_size = 1
            while True:
                alpha_t = alpha_g - step_size * jacobian[0]
                compressL_t = compressL_g - step_size * jacobian[1]

                suffix = f" -- mu {mu_t} -- totalMass {rho_t} -- angleRight {alpha_t} -- compressRatio {compressL_t}"
                cmd = cmd + suffix
                BASimProcess = subprocess.Popen(cmd, shell=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
                BASimProcess.communicate()
                fileName = f"datafiles/simDiscreteNet_l2_0.010000_compressRatio_{compressL_t:.6f}_angleRight_{alpha_t:.6f}_mu_{mu_t:.6f}_mass_{rho_t:.6f}.txt"
                data1 = np.loadtxt(fileName)
                idx = np.argmax(data1[:, 2])
                node= data1[idx]
                y_e = np.linalg.norm(target_p - node[1:])
                os.remove(fileName)
                if y_e < y_g or step_size < 1e-3:
                    break
                step_size *= 0.5

            alpha_g = alpha_g - step_size * jacobian[0]
            compressL_g = compressL_g - step_size * jacobian[1]
            y_g = y_e

            print(alpha_g, compressL_g, y_g, step_size)
            if step_size < 1e-3:
                break
            iter += 1
        print("==============")
        accuracy += y_g
        cache_result.append([alpha_g, compressL_g, y_g])

    elapsed_time = time.time() - start_time
    current, peak = tracemalloc.get_traced_memory()
    tracemalloc.stop()
    peak_memory = peak/(1024 ** 2)
    results = {"accuracy": accuracy/len(test_samples), "memory_usage": peak_memory, "time_cost": elapsed_time/len(test_samples)}
    print(results)

    return results
