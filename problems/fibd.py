import sys
import os

parent_directory = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.append(parent_directory)

from bioinphilmatics import core

print(core.mortal_fibonacci(82, 18))