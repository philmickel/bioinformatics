import sys
import os

parent_directory = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.append(parent_directory)

from bioinphilmatics import core

sample_dataset = [5, 3]
sample_output = 19

problem_dataset = [35, 5]
print(core.fibonacci(n=problem_dataset[0], k=problem_dataset[1]))