import matplotlib.pyplot as plt
from absl.testing import absltest
import numpy as np
import argparse
import sys
import subprocess
import os

def parse_custom_arguments():
    parser = argparse.ArgumentParser(description="Process custom key-value arguments.")
    parser.add_argument('--parameter', choices=['delta_alpha', 'epsilon', 'mu', 'mass'],
                        required=True, help="Specify the type of test to run.")
    args, remaining_args = parser.parse_known_args()


    return args, remaining_args

class TestOptimizationFunction(absltest.TestCase):
    def test_runs_without_error(self):
        global PARAMETER, testcases, base_node, results
        if PARAMETER == "delta_alpha":
            delta_alpha = testcases[0]
            epsilon = testcases[1]
            mu = testcases[2]
            bar_mass = testcases[3]
            EI = 200e9 * (0.05e-3)**3/ 12 * 5e-3
            mass = bar_mass * EI / 10 / 0.02**2
            delta_alpha_tests = np.linspace(delta_alpha * 0.8, delta_alpha * 1.2, 10)

            # compute  baseline first
            cmd = './simulations/simDER ./simulations/option.txt'
            suffix = f" -- mu {mu} -- totalMass {mass} -- angleRight {delta_alpha} -- compressRatio {1 - epsilon}"
            cmd = cmd + suffix
            BASimProcess = subprocess.Popen(cmd, shell=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
            BASimProcess.communicate()
            fileName = f"datafiles/simDiscreteNet_l2_0.010000_compressRatio_{1-epsilon:.6f}_angleRight_{delta_alpha:.6f}_mu_{mu:.6f}_mass_{mass:.6f}.txt"
            data1 = np.loadtxt(fileName)
            idx = np.argmax(data1[:, 2])
            base_node= data1[idx]
            base_node = base_node[1:]
            os.remove(fileName)

            print(base_node)
            results = {}
            for delta_alpha in delta_alpha_tests:
                cmd = './simulations/simDER ./simulations/option.txt'
                suffix = f" -- mu {mu} -- totalMass {mass} -- angleRight {delta_alpha} -- compressRatio {1 - epsilon}"
                cmd = cmd + suffix
                BASimProcess = subprocess.Popen(cmd, shell=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
                BASimProcess.communicate()
                fileName = f"datafiles/simDiscreteNet_l2_0.010000_compressRatio_{1-epsilon:.6f}_angleRight_{delta_alpha:.6f}_mu_{mu:.6f}_mass_{mass:.6f}.txt"
                data1 = np.loadtxt(fileName)
                idx = np.argmax(data1[:, 2])
                test_node= data1[idx]
                test_node = test_node[1:]
                results[delta_alpha] = test_node

        if PARAMETER == "epsilon":
            delta_alpha = testcases[0]
            epsilon = testcases[1]
            mu = testcases[2]
            bar_mass = testcases[3]
            epsilon_tests = np.linspace(epsilon * 0.8, epsilon * 1.2, 10)
            EI = 200e9 * (0.05e-3)**3/ 12 * 5e-3
            mass = bar_mass * EI / 10 / 0.02**2

            # compute  baseline first
            cmd = './simulations/simDER ./simulations/option.txt'
            suffix = f" -- mu {mu} -- totalMass {mass} -- angleRight {delta_alpha} -- compressRatio {1 - epsilon}"
            cmd = cmd + suffix
            BASimProcess = subprocess.Popen(cmd, shell=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
            BASimProcess.communicate()
            fileName = f"datafiles/simDiscreteNet_l2_0.010000_compressRatio_{1-epsilon:.6f}_angleRight_{delta_alpha:.6f}_mu_{mu:.6f}_mass_{mass:.6f}.txt"
            data1 = np.loadtxt(fileName)
            idx = np.argmax(data1[:, 2])
            base_node= data1[idx]
            base_node = base_node[1:]
            os.remove(fileName)

            results = {}
            for epsilon in epsilon_tests:
                cmd = './simulations/simDER ./simulations/option.txt'
                suffix = f" -- mu {mu} -- totalMass {mass} -- angleRight {delta_alpha} -- compressRatio {1 - epsilon}"
                cmd = cmd + suffix
                BASimProcess = subprocess.Popen(cmd, shell=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
                BASimProcess.communicate()
                fileName = f"datafiles/simDiscreteNet_l2_0.010000_compressRatio_{1-epsilon:.6f}_angleRight_{delta_alpha:.6f}_mu_{mu:.6f}_mass_{mass:.6f}.txt"
                data1 = np.loadtxt(fileName)
                idx = np.argmax(data1[:, 2])
                test_node= data1[idx]
                test_node = test_node[1:]
                results[epsilon] = test_node

        if PARAMETER == "mu":
            delta_alpha = testcases[0]
            epsilon = testcases[1]
            mu = testcases[2]
            bar_mass = testcases[3]
            mu_tests = np.linspace(mu * 0.8, mu * 1.2, 10)
            EI = 200e9 * (0.05e-3)**3/ 12 * 5e-3
            mass = bar_mass * EI / 10 / 0.02**2

            # compute  baseline first
            cmd = './simulations/simDER ./simulations/option.txt'
            suffix = f" -- mu {mu} -- totalMass {mass} -- angleRight {delta_alpha} -- compressRatio {1 - epsilon}"
            cmd = cmd + suffix
            BASimProcess = subprocess.Popen(cmd, shell=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
            BASimProcess.communicate()
            fileName = f"datafiles/simDiscreteNet_l2_0.010000_compressRatio_{1-epsilon:.6f}_angleRight_{delta_alpha:.6f}_mu_{mu:.6f}_mass_{mass:.6f}.txt"
            data1 = np.loadtxt(fileName)
            idx = np.argmax(data1[:, 2])
            base_node= data1[idx]
            base_node = base_node[1:]
            os.remove(fileName)

            results = {}
            for mu in mu_tests:
                cmd = './simulations/simDER ./simulations/option.txt'
                suffix = f" -- mu {mu} -- totalMass {mass} -- angleRight {delta_alpha} -- compressRatio {1 - epsilon}"
                cmd = cmd + suffix
                BASimProcess = subprocess.Popen(cmd, shell=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
                BASimProcess.communicate()
                fileName = f"datafiles/simDiscreteNet_l2_0.010000_compressRatio_{1-epsilon:.6f}_angleRight_{delta_alpha:.6f}_mu_{mu:.6f}_mass_{mass:.6f}.txt"
                data1 = np.loadtxt(fileName)
                idx = np.argmax(data1[:, 2])
                test_node= data1[idx]
                test_node = test_node[1:]
                results[mu] = test_node

        if PARAMETER == "mass":
            delta_alpha = testcases[0]
            epsilon = testcases[1]
            mu = testcases[2]
            bar_mass = testcases[3]
            bar_mass_tests = np.linspace(bar_mass * 0.8, bar_mass * 1.2, 10)

            EI = 200e9 * (0.05e-3)**3/ 12 * 5e-3
            mass = bar_mass * EI / 10 / 0.02**2

            # compute  baseline first
            cmd = './simulations/simDER ./simulations/option.txt'
            suffix = f" -- mu {mu} -- totalMass {mass} -- angleRight {delta_alpha} -- compressRatio {1 - epsilon}"
            cmd = cmd + suffix
            BASimProcess = subprocess.Popen(cmd, shell=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
            BASimProcess.communicate()
            fileName = f"datafiles/simDiscreteNet_l2_0.010000_compressRatio_{1-epsilon:.6f}_angleRight_{delta_alpha:.6f}_mu_{mu:.6f}_mass_{mass:.6f}.txt"
            data1 = np.loadtxt(fileName)
            idx = np.argmax(data1[:, 2])
            base_node= data1[idx]
            base_node = base_node[1:]
            os.remove(fileName)

            results = {}
            for bar_mass in bar_mass_tests:
                EI = 200e9 * (0.05e-3)**3/ 12 * 5e-3
                mass = bar_mass * EI / 10 / 0.02**2
                cmd = './simulations/simDER ./simulations/option.txt'
                suffix = f" -- mu {mu} -- totalMass {mass} -- angleRight {delta_alpha} -- compressRatio {1 - epsilon}"
                cmd = cmd + suffix
                BASimProcess = subprocess.Popen(cmd, shell=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
                BASimProcess.communicate()
                fileName = f"datafiles/simDiscreteNet_l2_0.010000_compressRatio_{1-epsilon:.6f}_angleRight_{delta_alpha:.6f}_mu_{mu:.6f}_mass_{mass:.6f}.txt"
                data1 = np.loadtxt(fileName)
                idx = np.argmax(data1[:, 2])
                test_node= data1[idx]
                test_node = test_node[1:]
                results[bar_mass] = test_node

        print(results)
        # process the data
        imperfection = np.linspace(-0.2, 0.2, 10)
        x_error = []
        y_error = []
        for key in results.keys():
            x_error.append(abs(results[key][0] - base_node[0])/base_node[0])
            y_error.append(abs(results[key][1] - base_node[1])/base_node[1])

        # create one figure with two subplots
        fig, ax = plt.subplots(2)
        ax[0].plot(imperfection, x_error, 'o-')
        ax[1].plot(imperfection, y_error, '^-')
        # add labels
        ax[0].set_xlabel('Imperfection ratio, %')
        ax[1].set_xlabel('Imperfection ratio, %')
        ax[0].set_ylabel('Realtive X error, %')
        ax[1].set_ylabel('Relative Y error, %')
        plt.show()

if __name__=="__main__":
    args, remaining_args = parse_custom_arguments()
    PARAMETER = args.parameter
    sys.argv = [sys.argv[0]] + remaining_args
    testcases = [0.1, 0.2, 0.3, 0.768] # delta_alpha, epsilon, mu, bar_mass

    absltest.main()
    print(base_node)
    print(results)
