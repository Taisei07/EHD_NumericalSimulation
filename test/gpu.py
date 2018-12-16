# -*- coding: utf-8 -*-

import cupy as cp
A = cp.arange(9).reshape(3, 3).astype('f') #cupy上で3*3の行列を生成
B = cp.arange(9).reshape(3, 3).astype('f') #cupy上で3*3の行列を生成
print('A = \n', A)
print('B = \n', B)
