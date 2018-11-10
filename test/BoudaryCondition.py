# -*- coding: utf-8 -*-

import numpy as np
import matplotlib.pyplot as plt

#境界条件を設定する

#格子定義点
u = np.array([[0] * ms for i in range(n)])#ms:x方向,n:y方向
v = np.array([[0] * ms for i in range(n)])
p = np.array([[0] * ms for i in range(n)])
