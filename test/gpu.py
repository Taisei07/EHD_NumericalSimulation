# -*- coding: utf-8 -*-

import numpy as np
import cupy as cp

x_cpu = np.zeros((10,10))

x_gpu = cp.asarray(x_cpu)
x_cpu = cp.asnumpy(x_gpu)

print(x_gpu ** 2)

xp = cp.get_array_module(x_gpu)
print(xp.square(x_gpu))
