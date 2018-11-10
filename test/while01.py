# -*- coding: utf-8 -*-

import numpy as np
import matplotlib.pyplot as plt

#初期値
deltaT = 1
deltax = 1
deltay = 1
#格子数
ms = 6#x方向は600格子。座標番号は0~599
n = 7#y方向は200格子。座標番号は0~199
#格子定義点
u = np.array([[0] * ms for i in range(n)])#ms:x方向,n:y方向
v = np.array([[0] * ms for i in range(n)])
p = np.array([[0] * ms for i in range(n)])
#物性値
rho = 1520 #HFE-7100
nu = 0.38 #HFE-7100
#開始格子点の指定
i = 0
j = 0
#計算開始
while 0 <= j <= (n-2):
    while 0 <= i <= (ms-2):
        u[i][j] = 2
        i += 1
    j += 1
    i = 1
print u
