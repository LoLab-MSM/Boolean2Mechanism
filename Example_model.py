
# Boolean T-cell model from Saez-Rodriguez et al., 2007

from boolean2 import Model
import numpy as np
import matplotlib.pyplot as plt

rules = """
A = False
B = False
C = True

A* = B and C
B* = B
C* = C
"""
