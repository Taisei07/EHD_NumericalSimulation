# -*- coding: utf-8 -*-

import numpy as np
import matplotlib.pyplot as plt

#deltapの設定用

#物性値
rho = 1000 #密度[kg/m^3]
nu = 1.004#動粘度[m^2/s]
#定数
ms = 60#x方向の要素数
n = 20#y方向の要素数
deltaT = 1#時間刻み
deltax = 1#x方向の要素間距離
deltay = 1#y方向の要素間距離
omega = 0.5#圧力補正における緩和係数
M = 1#連続の式の収束条件
#境界条件
p_boudaryAH = 100#入り口圧力[Pa]
u_boundaryAH = 30#入り口速度[m/s]

#格子定義点(oldはtimestep=n, newはtimestep=n+1)
u_old = list([[0] * n for i in range(ms)])#ms:x方向,n:y方向
v_old = list([[0] * n for i in range(ms)])
p_old = list([[0] * n for i in range(ms)])
u_new = list([[0] * n for i in range(ms)])
v_new = list([[0] * n for i in range(ms)])
p_new = list([[0] * n for i in range(ms)])

#i,jの初期化
i = 0
j = 1

#境界条件のループ
while 0 < j < n-1:
    p_old[i][j] = p_boudaryAH
    u_old[i][j] = u_boundaryAH
    j += 1
print "boudary condition setting"

#i,jの初期化
i = 1
j = 1

#deltap,deltau,deltavの計算
while 0 < i < ms-1:
    while 0 < j < n-1:
        deltap = -rho / (2*deltaT) * (deltax * deltay) / (deltax**2 + deltay**2) * (deltay * (u_old[i][j] - u_old[i-1][j]) + deltax * (v_old[i][j] - v_old[i][j-1]))
        print deltap
        p_old[i][j] = p_old[i][j] + (omega * deltap)
        u_old[i][j] = u_old[i][j] + 1 / rho * deltaT / deltax * deltap
        v_old[i][j] = v_old[i][j] + 1 / rho * deltaT / deltay * deltap
        j += 1
    j = 1
    i += 1
print "end"
print p_old
