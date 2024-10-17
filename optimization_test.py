from absl import flags
from absl.testing import absltest
from absl.testing import parameterized
import numpy as np
import tracemalloc
import argparse
import sys
import utils

def parse_custom_arguments():
    parser = argparse.ArgumentParser(description="Process custom key-value arguments.")
    parser.add_argument('--test_type', choices=['data-driven', 'bayesian-opt', 'gradient-descent'],
                        required=True, help="Specify the type of test to run.")
    args, remaining_args = parser.parse_known_args()


    return args, remaining_args


class TestOptimizationFunction(absltest.TestCase):
    def test_runs_without_error(self):
        global TEST_TYPE, testcases
        if TEST_TYPE == "data-driven":
            results = utils.data_driven(testcases)

        if TEST_TYPE == "bayesian-opt":
            results = utils.bayestian_opt(testcases)

        if TEST_TYPE == "gradient-descent":
            results = utils.gradient_descent(testcases)

        print(results)

def load_testcases():
    data = np.loadtxt("test_cases.txt").tolist()
    return data


if __name__=="__main__":
    args, remaining_args = parse_custom_arguments()
    TEST_TYPE = args.test_type
    sys.argv = [sys.argv[0]] + remaining_args
    testcases = load_testcases()

    absltest.main()

    print(TEST_TYPE)
    
