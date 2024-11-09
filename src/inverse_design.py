import numpy as np
import matplotlib.pyplot as plt
from utils import loadModel, loadData, getRange, pred, inverseOptimize,  checkTarget
import subprocess

import argparse
import time
import os
import tracemalloc

def parse_custom_arguments():
    parser = argparse.ArgumentParser(description="Process custom key-value arguments.")
    parser.add_argument("args", nargs='*', help="key=value pairs")

    args = parser.parse_args()
    arg_dict = {}
    
    for arg in args.args:
        if ":=" in arg:
            key, value = arg.split(":=", 1)
            if value.isdigit():
                value = int(value)
            elif value.lower() == 'true':
                value = True
            elif value.lower() == 'false':
                value = False
            arg_dict[key] = value
        else:
            raise ValueError(f"Argument {arg} is not in the key:=value format.")
    
    return arg_dict
    
def main():
    arguments = parse_custom_arguments()
    print("Parsed arguments:", arguments)
    test_num = arguments.get('test_num')
    plot_flag = arguments.get('plot')

    # load model
    yaml_path = "model.yaml"
    checkpoint_path = "./output/checkpoints/checkpoint_200.pt"
    model = loadModel(yaml_path=yaml_path, checkpoint_path=checkpoint_path)

    # load data
    x, y, mu, rho, alpha, compressL = loadData()
    alpha_f, alpha_c = getRange(alpha)
    compressL_f, compressL_c = getRange(compressL)
    rho_f, rho_c = getRange(rho)
    mu_f, mu_c = getRange(mu)

    Data = []

    iter = 0
    while iter < test_num:
        # determine the parameters
        mu_rand = np.random.uniform(mu_f, mu_c)
        rho_rand = np.random.uniform(rho_f, rho_c)
        # check the available region
        alpha = np.linspace(alpha_f, alpha_c, 10)
        compressL = np.linspace(compressL_f, compressL_c, 10)

        alpha, compressL= np.meshgrid(alpha, compressL, indexing='ij')
        alpha_input = alpha.flatten()
        compressL_input = compressL.flatten()
        mu_input = np.ones_like(alpha_input) * mu_rand
        rho_input = np.ones_like(alpha_input) * rho_rand

        features = np.vstack([alpha_input, compressL_input, mu_input, rho_input]).T
        outputs = pred(model, features)


        # generate a random point in the region
        x_rand = np.random.uniform(np.min(outputs[:, 0]), np.max(outputs[:, 0]))
        y_rand = np.random.uniform(np.min(outputs[:, 1]), np.max(outputs[:, 1]))
        target_p = np.asarray([x_rand, y_rand]).reshape(-1, 2)

        flag = checkTarget(outputs, target_p, plot_flag)

        if not flag:
            print("The target is not achievable")
            if plot_flag:
                plt.plot(target_p[0, 0], target_p[0, 1], 'x',  label = "Target")
                plt.legend()
                plt.title(rf"$\mu$ = {mu_rand:.2g}, $\rho$= {rho_rand:.2g}")
                plt.show()
            continue

        # selecet the initial guess
        offset = outputs - target_p
        offset = np.sqrt((offset**2).sum(axis = 1))

        min_index = np.argmin(offset)

        alpha_int = alpha_input[min_index]
        compressL_int = compressL_input[min_index]

        inputs = [alpha_int, compressL_int, mu_rand, rho_rand]
        alpha, compressL, mu, rho = inverseOptimize(inputs, model, target_p.copy())

        cmd = './simulations/simDER ./simulations/option.txt'
        suffix = f" -- mu {mu} -- totalMass {rho} -- angleRight {alpha} -- compressRatio {compressL}"
        cmd = cmd + suffix
        BASimProcess = subprocess.Popen(cmd, shell=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
        BASimProcess.communicate()
        fileName = f"datafiles/simDiscreteNet_l2_0.010000_compressRatio_{compressL:.6f}_angleRight_{alpha:.6f}_mu_{mu:.6f}_mass_{rho:.6f}.txt"
        data1 = np.loadtxt(fileName)
        idx = np.argmax(data1[:, 2])
        node= data1[idx]
        offset = np.linalg.norm(target_p - node[1:])

        Data.append(offset)

        iter += 1
        if plot_flag:
            plt.plot(data1[:, 1], data1[:, 2], label = 'Path')
            plt.plot(target_p[0, 0], target_p[0, 1], 'x',  label = "Target")
            plt.legend()
            plt.title(rf"$\mu$ = {mu_rand:.2g}, $\rho$= {rho_rand:.2g}, $\delta \alpha$ = {alpha:.2g}, $\eta$ = {1-compressL:.2g}")
            plt.axis('equal')
            plt.show()
        
        os.remove(fileName)

    Data = np.asarray(Data)
    print(f"The average error of the inverse design is {np.mean(Data)}, and standard deviation is {np.std(Data)}")

if __name__ == "__main__":
    main()
