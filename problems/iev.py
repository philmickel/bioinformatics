import sys
import os
import numpy as np

parent_directory = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.append(parent_directory)

from bioinphilmatics import core

sample_dataset = "1 0 0 1 0 1"
sample_output = 3.5

def problem_formatting(dataset):
    # AA-AA AA-Aa AA-aa Aa-Aa Aa-aa aa-aa
    PROBABILITIES = np.array([1, 1, 1, 0.75, 0.5, 0])
    dataset_list = dataset.split()
    int_list = list(map(int, dataset_list))
    dataset_array = np.array(int_list)
    return 2*np.dot(PROBABILITIES, dataset_array)

problem_dataset = "17642 16582 16324 16281 17686 17301"

print(problem_formatting(problem_dataset))