# -*- coding: utf-8 -*-
"""
Created on Thu Apr 11 11:09:33 2024

@author: cyril
"""

import numpy as np
import matplotlib.pyplot as plt

m = np.random.rand(100, 100)
print(m)

for row in m:
    pos = np.where(row < 0.01)
    print(m[:, pos[0]])