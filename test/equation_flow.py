# -*- coding: utf-8 -*-

import numpy as np
import matplotlib.pyplot as plt

#初期値
deltaT = 1
deltax = 1
deltay = 1
#格子数
ms = 100
n = 100
#格子定義点
u = [[0] * ms for i in range(n)]
v = [[0] * ms for i in range(n)]
p = [[0] * ms for i in range(n)]
#物性値
rho = 1520 #HFE-7100
nu = 0.38 #HFE-7100
#計算格子の始まり
i = 0
j = 0
#計算開始
while i < ms:
    while j < n:
        u[i][j] = u[i][j] + \
                deltaT*(-(u[i][j]*(u[i+1][j]-u[i-1][j])/(2*deltax) - \
                (abs(u[i][j])*(u[i+1][j]-2*u[i][j]+u[i-1][j])/(2*deltax)) - \
                (v[i][j]*(u[i][j+1]-u[i][j-1])/(2*deltay)-abs(v[i][j])*(u[i][j+1]-2*u[i][j]+u[i][j-1])/(2*deltay)) - \
                ((1/rho)*(p[i+1][j])/(deltax)) + \
                nu((u[i+1][j] - 2*u[i][j] + u[i][j-1]) / ((deltax)**2) + (u[i][j+1] - 2*u[i][j] + u[i][j-1]) / ((deltay)**2))
        j += 1
    print "j = " + str(j)
    i += 1
    j = 0
print "end"
