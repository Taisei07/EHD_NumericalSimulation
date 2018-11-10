# -*- coding: utf-8 -*-

import numpy as np
import matplotlib.pyplot as plt

#圧力補正量を設定する

#物性値
rho = 950.0 #密度[kg/m^3]
nu = 1.004#動粘度[m^2/s]
#定数
ms = 60#x方向の要素数
n = 20#y方向の要素数
deltaT = 1.0#時間刻み
deltax = 1.0#x方向の要素間距離
deltay = 1.0#y方向の要素間距離
omega = 0.5#圧力補正における緩和係数
M = 1#連続の式の収束条件
#境界条件
u_boundaryAH = 0.0#入り口x方向速度[m/s]
u_boundaryFG = 0.0#出口x方向速度[m/s]
p_boundaryAH = 500.0#入り口圧力[Pa]
p_boundaryFG = 0.0#出口圧力[Pa]
#格子定義点(oldはtimestep=n, newはtimestep=n+1)
u_old = np.array([[0.0] * n for i in range(ms)])#ms:x方向,n:y方向
v_old = np.array([[0.0] * n for i in range(ms)])
p_old = np.array([[0.0] * n for i in range(ms)])
u_new = np.array([[0.0] * n for i in range(ms)])
v_new = np.array([[0.0] * n for i in range(ms)])
p_new = np.array([[0.0] * n for i in range(ms)])

#境界条件の設定
#入口境界
i = 0
j = 0
while j <= n-1:
    u_old[0][j] = u_boundaryAH
    u_old[1][j] = u_boundaryAH
    p_old[0][j] = p_boundaryAH
    p_old[1][j] = p_boundaryAH
    j += 1
print "入り口境界設定完了"
#出口境界
i = 0
j = 0
while j <= n-1:
    u_old[ms-1][j] = u_boundaryFG
    u_old[ms-2][j] = u_boundaryFG
    p_old[ms-1][j] = p_boundaryFG
    p_old[ms-2][j] = p_boundaryFG
    j += 1
print "出口境界設定完了"

#境界条件の確認
print u_old
print v_old
print p_old

print "ここから計算開始"

#連続の式の定義
#DIV = (u_old[i][j] - u_old[i-1][j])/deltax + (v_old[i][j] - v_old[i][j-1])/deltay
#対流項の定義
#CNVU = (u_old[i][j] * (u_old[i+1][j]-u_old[i-1][j]) /(2*deltax) - abs(u_old[i][j]) * (u_old[i+1][j]-2*u_old[i][j]+u_old[i-1][j]) / (2*deltax)) - (v_old[i][j]*(u_old[i][j+1] - u_old[i][j-1]) / (2 * deltay) - abs(v_old[i][j])*(u_old[i][j+1]-2*u_old[i][j]+u_old[i][j-1]) / deltay)
#CNVV = (u_old[i][j] * (v_old[i+1][j]-v_old[i-1][j]) /(2*deltax) - abs(u_old[i][j]) * (v_old[i+1][j]-2*v_old[i][j]+v_old[i-1][j]) / (2*deltax)) - (v_old[i][j]*(v_old[i][j+1] - v_old[i][j-1]) / (2 * deltay) - abs(v_old[i][j])*(v_old[i][j+1]-2*v_old[i][j]+v_old[i][j-1]) / deltay)
#粘性項
#DIFU = nu * ((u_old[i+1][j]-2*u_old[i][j]+u_old[i-1][j]) / (deltax)**2 + (u_old[i][j+1]-2*u_old[i][j]+u_old[i][j-1]) / (deltay)**2)
#DIFV = nu * ((v_old[i+1][j]-2*v_old[i][j]+v_old[i-1][j]) / (deltax)**2 + (v_old[i][j+1]-2*v_old[i][j]+v_old[i][j-1]) / (deltay)**2)
#圧力修正量deltapと速度場修正量の定義
#deltap = -rho/(2*deltaT)*(deltax * deltay)/(deltax**2 + deltay**2)*(deltay*(u_old[i][j]-u_old[i-1][j])+deltax*(v_old[i][j]-v_old[i][j-1]))
#deltau = 1/rho*deltaT/deltax*deltap
#deltav = 1/rho*deltaT/deltay*deltap

#u_old,v_oldの仮値を設定
i = 1
j = 1
while 0 < i < ms-1:
    while 0 < j < n-1:
        CNVU = (u_old[i][j] * (u_old[i+1][j]-u_old[i-1][j]) * 1.0 /(2*deltax) - abs(u_old[i][j]) * (u_old[i+1][j]-2*u_old[i][j]+u_old[i-1][j]) * 1.0 / (2*deltax)) - (v_old[i][j]*(u_old[i][j+1] - u_old[i][j-1]) * 1.0 / (2 * deltay) - abs(v_old[i][j])*(u_old[i][j+1]-2*u_old[i][j]+u_old[i][j-1]) * 1.0 / deltay)
        CNVV = (u_old[i][j] * (v_old[i+1][j]-v_old[i-1][j]) * 1.0 /(2*deltax) - abs(u_old[i][j]) * (v_old[i+1][j]-2*v_old[i][j]+v_old[i-1][j]) * 1.0 / (2*deltax)) - (v_old[i][j]*(v_old[i][j+1] - v_old[i][j-1]) * 1.0 / (2 * deltay) - abs(v_old[i][j])*(v_old[i][j+1]-2*v_old[i][j]+v_old[i][j-1]) * 1.0 / deltay)
        DIFU = nu * ((u_old[i+1][j]-2*u_old[i][j]+u_old[i-1][j]) * 1.0 / ((deltax)**2) + (u_old[i][j+1]-2*u_old[i][j]+u_old[i][j-1]) * 1.0 / ((deltay)**2))
        DIFV = nu * ((v_old[i+1][j]-2*v_old[i][j]+v_old[i-1][j]) * 1.0 / ((deltax)**2) + (v_old[i][j+1]-2*v_old[i][j]+v_old[i][j-1]) * 1.0 / ((deltay)**2))
        u_old[i][j] = deltaT * (u_old[i][j] - (1.0/rho)*(p_old[i+1][j]-p_old[i][j])/deltax - CNVU + DIFU)
        v_old[i][j] = deltaT * (v_old[i][j] - (1.0/rho)*(p_old[i][j+1]-p_old[i][j])/deltay - CNVV + DIFV)
        print "y方向にループ計算中"
        print "u_old[" + str(i) + "]" + "[" + str(j) + "] = " + str(u_old[i][j])
        print "v_old[" + str(i) + "]" + "[" + str(j) + "] = " + str(v_old[i][j])
        print "CNVU = " + str(CNVU)
        print "CNVV = " + str(CNVV)
        print "DIFU = " + str(DIFU)
        print "DIFV = " + str(DIFV)
        j += 1
    j = 1
    i += 1
    print "x方向にループ計算中"
print u_old
print v_old
print "u,vの仮値設定完了"
