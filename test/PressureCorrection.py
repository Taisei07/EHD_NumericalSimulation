# -*- coding: utf-8 -*-

import numpy as np
import matplotlib.pyplot as plt

#圧力補正量を設定する

#物性値
rho = 1000 #[kg/m^3]
nu = 0.01
#定数
ms = 10
n = 10
deltaT = 0.001
deltax = 1
deltay = 1
omega = 0.5
#境界条件
u_boudaryAH = 10

#格子定義点(oldはtimestep=n, newはtimestep=n+1)
u_old = np.array([[0] * ms for i in range(n)])#ms:x方向,n:y方向
v_old = np.array([[0] * ms for i in range(n)])
p_old = np.array([[0] * ms for i in range(n)])
u_new = np.array([[0] * ms for i in range(n)])
v_new = np.array([[0] * ms for i in range(n)])
p_new = np.array([[0] * ms for i in range(n)])

#iとjをリセット
i = 0
j = 0

#境界条件の設定
#入口速度を設定する
while j <= n-1:
    u_old[0][j] = u_boudaryAH
    j += 1
#境界条件の確認
print u_old
print v_old
print p_old

print "ここから計算開始"

#iとjをリセット
i = 0
j = 0

#連続の式の定義
DIV = (u_new[i][j] - u_new[i-1][j])/deltax + (v_new[i][j] - v_new[i][j-1])/deltay
#対流項の定義
CNVU = (u_old[i][j] * (u_old[i+1][j]-u_old[i-1][j]) /(2*deltax) - abs(u_old[i][j]) * (u_old[i+1][j]-2*u_old[i][j]+u_old[i-1][j]) / (2*deltax)) - (v_old[i][j]*(u_old[i][j+1] - u_old[i][j-1]) / (2 * deltay) - abs(v_old[i][j])*(u_old[i][j+1]-2*u_old[i][j]+u_old[i][j-1]) / deltay)
CNVV = (u_old[i][j] * (v_old[i+1][j]-v_old[i-1][j]) /(2*deltax) - abs(u_old[i][j]) * (v_old[i+1][j]-2*v_old[i][j]+v_old[i-1][j]) / (2*deltax)) - (v_old[i][j]*(v_old[i][j+1] - v_old[i][j-1]) / (2 * deltay) - abs(v_old[i][j])*(v_old[i][j+1]-2*v_old[i][j]+v_old[i][j-1]) / deltay)
#粘性項
DIFU = nu * ((u_old[i+1][j]-2*u_old[i][j]+u_old[i-1][j]) / (deltax)**2 + (u_old[i][j+1]-2*u_old[i][j]+u_old[i][j-1]) / (deltay)**2)
DIFV = nu * ((v_old[i+1][j]-2*v_old[i][j]+v_old[i-1][j]) / (deltax)**2 + (v_old[i][j+1]-2*v_old[i][j]+v_old[i][j-1]) / (deltay)**2)
#圧力修正量deltapと速度場修正量の定義
deltap = -rho/(2*deltaT)*(deltax * deltay)/(deltax**2 + deltay**2)*(deltay*(u_old[i][j]-u_old[i-1][j])+deltax*(v_old[i][j]-v_old[i][j-1]))
deltau = 1/rho*deltaT/deltax*deltap
deltav = 1/rho*deltaT/deltay*deltap

#iとjの値をリセット

#u_old,v_oldの仮値を設定
while 0 < j < n-1:
    while 0 < i < ms-1:
        CNVU = (u_old[i][j] * (u_old[i+1][j]-u_old[i-1][j]) /(2*deltax) - abs(u_old[i][j]) * (u_old[i+1][j]-2*u_old[i][j]+u_old[i-1][j]) / (2*deltax)) - (v_old[i][j]*(u_old[i][j+1] - u_old[i][j-1]) / (2 * deltay) - abs(v_old[i][j])*(u_old[i][j+1]-2*u_old[i][j]+u_old[i][j-1]) / deltay)
        CNVV = (u_old[i][j] * (v_old[i+1][j]-v_old[i-1][j]) /(2*deltax) - abs(u_old[i][j]) * (v_old[i+1][j]-2*v_old[i][j]+v_old[i-1][j]) / (2*deltax)) - (v_old[i][j]*(v_old[i][j+1] - v_old[i][j-1]) / (2 * deltay) - abs(v_old[i][j])*(v_old[i][j+1]-2*v_old[i][j]+v_old[i][j-1]) / deltay)
        DIFU = nu * ((u_old[i+1][j]-2*u_old[i][j]+u_old[i-1][j]) / (deltax)**2 + (u_old[i][j+1]-2*u_old[i][j]+u_old[i][j-1]) / (deltay)**2)
        DIFV = nu * ((v_old[i+1][j]-2*v_old[i][j]+v_old[i-1][j]) / (deltax)**2 + (v_old[i][j+1]-2*v_old[i][j]+v_old[i][j-1]) / (deltay)**2)
        u_old[i][j] = deltaT * (u_old[i][j] - (1/rho)*(p_old[i+1][j]-p_old[i][j])/deltax - CNVU + DIFU)
        v_old[i][j] = deltaT * (v_old[i][j] - (1/rho)*(p_old[i][j+1]-p_old[i][j])/deltay - CNVV + DIFV)
        i += 1
    i = 1
    j += 1
print u_old
print v_old
#圧力補正のループ
while DIV >  0.1:#連続の式の収束条件を満たさない場合
    u_new[i][j] = u_old[i][j] + omega * deltau
    v_new[i][j] = v_old[i][j] + omega * deltav
    p_new[i][j] = p_old[i][j] + omega * deltap
    u_old[i][j] = u_new[i][j]
    v_old[i][j] = v_new[i][j]
    p_old[i][j] = p_new[i][j]
    print "DIV = " + str(DIV)
print "u_new[1][1] = " + str(u_new[1][1])

#各変数の確認
print "DIV = " + str(DIV)
print "CNVU = " + str(CNVU)
print "CNVV = " + str(CNVV)
print "DIFU = " + str(DIFU)
print "DIFV = " + str(DIFV)
print "deltap = " + str(deltap)
print "deltau = " + str(deltau)
print "deltav = " + str(deltav)
#速度場、圧力場の確認
print u_new
print v_new
print p_new
