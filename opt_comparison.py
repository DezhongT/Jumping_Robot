import numpy as np
import matplotlib.pyplot as plt
from utils import loadModel, loadData, getRange, pred, inverseOptimize,  checkTarget
import subprocess

import argparse
import time

from argparse import Namespace
from dragonfly import load_config
from dragonfly.exd.experiment_caller import CPFunctionCaller
from dragonfly.opt import gp_bandit
import os


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
    test_num = 1
    plot_flag = True
    

    # inverse design: alpha 0.031080069977683103, compressL 0.7165367518034843, mu 0.4231008938838332 rho 0.001292034620946721


    # load model
    yaml_path = "model.yaml"
    checkpoint_path = "./output/checkpoints/forward_model.pt"
    model = loadModel(yaml_path=yaml_path, checkpoint_path=checkpoint_path)

    # load data
    x, y, mu, rho, alpha, compressL = loadData()
    alpha_f, alpha_c = getRange(alpha)
    compressL_f, compressL_c = getRange(compressL)
    rho_f, rho_c = getRange(rho)
    mu_f, mu_c = getRange(mu)

    print(alpha_c, alpha_f, compressL_f, compressL_c, mu_f, mu_c, rho_f, rho_c)


    mu =  0.4231008938838332
    rho = 0.001292034620946721
    target_p = np.array([0.15070193, 0.81060987]).reshape(-1, 2)

    start_time = time.time()
    # let us try bayesian optimization
    num_init = 6
    batch_size = 5
    
    options = Namespace(
        build_new_model_every = batch_size,
        init_captial = num_init-1,
        gpb_hp_tune_criterion = "ml-post_sampling",
    )
    domain_vars = [{'type': 'float',  'min': 0.01, 'max': 0.19, 'dim': 1}, # alpha
                   {'type': 'float',  'min': 0.7, 'max': 0.9, 'dim': 1}, # compressL
                  ]
    
    config_params = {"domain": domain_vars}
    config = load_config(config_params)
    func_caller = CPFunctionCaller(None, config.domain, domain_orderings=config.domain_orderings)
    opt = gp_bandit.CPGPBandit(func_caller, 'default', ask_tell_mode=True, options=options)  # opt is the optimizer object
    opt.initialise()

    cache = {}
    # warm up
    initial_samples = opt.ask(num_init)
    for sample in initial_samples:
        alpha = sample[0][0]
        compressL = sample[1][0]
        # check if the file exist
        if (alpha, compressL) in cache:
            opt.tell([(sample, -cache[(alpha, compressL)])])
        else:
            cmd = './simulations/simDER ./simulations/option.txt'
            suffix = f" -- mu {mu} -- totalMass {rho} -- angleRight {alpha} -- compressRatio {compressL}"
            cmd = cmd + suffix
            BASimProcess = subprocess.Popen(cmd, shell=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
            BASimProcess.communicate()
            fileName = f"datafiles/simDiscreteNet_l2_0.010000_compressRatio_{compressL:.6f}_angleRight_{alpha:.6f}_mu_{mu:.6f}_mass_{rho:.6f}.txt"
            data1 = np.loadtxt(fileName)
            idx = np.argmax(data1[:, 2])
            node= data1[idx]
            y = np.linalg.norm(target_p - node[1:])
            opt.tell([(sample, -y)])
            cache[(alpha, compressL)] = y
            os.remove(fileName)
            print(y)

    # update model
    opt._build_new_model()
    opt._set_next_gp()

    epochs = 0
    while epochs < 20:
        # sampling
        batch_samples = []
        for i in range(batch_size):
            batch_samples.append(opt.ask())

        average_err = 0
        for sample in batch_samples:
            alpha = sample[0][0]
            compressL = sample[1][0]
            # check if the file exist
            if (alpha, compressL) in cache:
                opt.tell([(sample, -cache[(alpha, compressL)])])
            else:
                cmd = './simulations/simDER ./simulations/option.txt'
                suffix = f" -- mu {mu} -- totalMass {rho} -- angleRight {alpha} -- compressRatio {compressL}"
                cmd = cmd + suffix
                BASimProcess = subprocess.Popen(cmd, shell=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
                BASimProcess.communicate()
                fileName = f"datafiles/simDiscreteNet_l2_0.010000_compressRatio_{compressL:.6f}_angleRight_{alpha:.6f}_mu_{mu:.6f}_mass_{rho:.6f}.txt"
                data1 = np.loadtxt(fileName)
                idx = np.argmax(data1[:, 2])
                node= data1[idx]
                y = np.linalg.norm(target_p - node[1:])
                opt.tell([(sample, -y)])
                cache[(alpha, compressL)] = y
                os.remove(fileName)
            average_err += cache[(alpha, compressL)]

        print(f"epoch {epochs}, error: {average_err/batch_size}")
        opt._build_new_model()
        opt._set_next_gp()
        epochs += 1
    print(f"Time elasped: {time.time() - start_time}")
    ## begin optimization
    print(cache)



    exit(0)

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

        print(f"inverse design: alpha {alpha}, compressL {compressL}, mu {mu} rho {rho}")
        print(target_p)

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

    Data = np.asarray(Data)
    print(f"The average error of the inverse design is {np.mean(Data)}, and standard deviation is {np.std(Data)}")

if __name__ == "__main__":
    main()
