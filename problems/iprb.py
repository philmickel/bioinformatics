import sys
import os

parent_directory = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.append(parent_directory)

from bioinphilmatics import core

sample_dataset = "2 2 2"
sample_output = 0.78333

problem_dataset = "26 27 23"

dataset_split = problem_dataset.split(' ')
print(core.mendellian_inheritance(k=int(dataset_split[0]), m=int(dataset_split[1]), n=int(dataset_split[2])))