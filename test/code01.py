# -*- coding: utf-8 -*-

import numpy as np
import matplotlib.pyplot as plt

#物性値
#粘性係数
nu = 0
#流体密度
rho = 0
#誘電率
epsilon = 0
#イオン移動度
K = 0
#導電率
sigma = 0
#イオン拡散係数
Di = 0
#初期値・境界条件の設定
#x方向の速度
uw = #BoundaryAHのx方向速度
un = #WallGHのx方向速度
ue = #BoundaryFGのx方向速度
us = #WallAFのx方向速度
#y方向の速度
vw = #BoundaryAHのy方向速度
ve = #BoundaryFGのy方向速度
#圧力
pw = #BoundaryAHの圧力
pn = #WallGHの圧力
pe = #BoundaryFGの圧力
ps = #WallAFの圧力
#電荷密度
qw = #BoundaryAHの電荷密度
qn = #WallGHの電荷密度
qe = #BoundaryFGの電荷密度
qs = #WallAFの電荷密度
#電位
phiw = #BoundaryAHの電位
phin = #WallGHの電位
phie = #BoundaryFGの電位
phis = #WallAFの電位
phiE1 = #ElectrodeBCの電位
phiE2 = #ElectrodeDEの電位
#時間の定義
deltaT = 1 #解析の刻み幅
T = 1
Tend = 20 #解析時間設定
#ステップ数
loop = 500
#流れ場の定義
i = 0 #i方向のメッシュ始まり
j = 0 #j方向のメッシュ始まり
ms = 599 #=(x方向メッシュ数-1)
n = 199 #=(y方向メッシュ数-1)
m1 = 189 #電極上B点の座標
m2 = 289 #電極上C点の座標
m3 = 309 #電極上D点の座標
m4 = 409 #電極上E点の座標
#格子のサイズ
deltax = 1 #x方向グリッドのサイズ
deltay = 1 #y方向グリッドのサイズ
#反復法での許容誤差
eps =
#緩和係数
beta =

#未知数の設定
#x方向速度（x方向ms、y方向nの座標をリストで設定完了。u[ms(x)][n(y)]で座標指定可能。）
u = [[0] * ms for i in range(n)]
#y方向の速度
v = [[0] * ms for i in range(n)]
#圧力
p = [[0] * ms for i in range(n)]
#電荷密度
q = [[0] * ms for i in range(n)]
#電位
phi = [[0] * ms for i in range(n)]

#境界条件の設定
#BoundaryAH
while j <= n:
    u[0][j] = uw
    u[1][j] = uw
    v[0][j] = vw
    v[1][j] = vw
    p[0][j] = pw
    p[1][j] = pw
    phi[0][j] = phiw
    phi[1][j] = phiw
    q[0][j] = qw
    q[1][j] = qw
    j = j + 1
print "u[0][j] = " + str(uw) + ", " + "v[0][j] = " + str(vw) + ", " +"p[0][j] = " + str(pw) + ", " + "phi[0][j] = " + str(phiw) + ", " + "q[0][j] = " + str(qw)

#WallGH
while i <= ms:
    u[i][n] = un
    u[i][n-1] = un
    v[i][n] = 0
    p[i][n] = pn
    p[i][n-1] = pn
    phi[i][n] = phin
    phi[i][n-1] = phin
    q[i][n] = qn
    q[i][n-1] = qn
    i = i + 1
print "u[i][n] = " + str(un) + ", " + "v[i][n] = 0" + ", " +"p[i][n] = " + str(pn) + ", " + "phi[i][n] = " + str(phin) + ", " + "q[i][n] = " + str(qn)

#BoundaryFG
j = 0
while j <= n:
    u[ms][j] = ue
    u[ms-1][j] = ue
    v[ms][j] = ve
    v[ms-1][j] = ve
    p[ms][j] = pe
    p[ms-1][j] = pe
    phi[ms][j] = phie
    phi[ms-1][j] = phie
    q[ms][j] = qe
    q[ms-1][j] = qe
    j = j + 1
print "u[ms][j] = " + str(ue) + ", " + "v[ms][j] = " + ", " + str(ve) + ", " + "p[ms][j] = " + str(pe) + ", " + "phi[ms][j] = " + str(phie) + ", " + "q[ms][j] = " + str(qe)

#WallAB,WallCD,WallEF（ElectrodeBC,DEの部分も設定してしまう）
i = 0
while i < ms:
    u[i][0] = -us
    u[i][1] = -us
    v[i][0] = 0
    p[i][0] = ps
    p[i][1] = ps
    phi[i][0] = phis
    phi[i][1] = phis
    q[i][0] = qs
    q[i][1] = qs
    i = i + 1
print "u[i][0] = " + str(-us) + ", " + "v[i][0] = 0 " + "p[i][0] = " + str(ps) + ", " + "phi[i][0] = " + str(phis) + ", " + "q[i][0] = " + str(qs)

#ElectrodeBC,DE
while m1 <= i < m2 :
    phi[i][0] = phiE1
    q[i][0] = -epsilon*(phi[i][2]-phi[i][1])/deltay**2
    i = i + 1
print "phi[i][0] = " + str(phiE1)
while m3 <= i < m4 :
    phi[i][0] = phiE2
    q[i][0] = -epsilon*(phi[i][2]-phi[i][1])/deltay**2
    i = i + 1
print "phi[i][0] = " + str(phiE2)

#計算開始
i = 1#境界条件で0は決まっているから1からスタートする？
j = 1
while T < Tend: #解析時間の指定
    while u2 + v2 < 10: #条件式(連続の式の条件を設定する)
        #uが正負の場合それぞれに対して計算処理をする=>abs()を使って絶対値を計算している
        u[i][j] = u[i][j] + ¥
                deltaT*(-(u[i][j] * (u[i+1][j] - u[i-1][j]) / (2*deltax) - ¥
                (abs(u[i][j]) * (u[i+1][j] - 2*u[i][j] + u[i-1][j]) / (2*deltax)) - ¥
                (v[i][j] * (u[i][j+1] - u[i][j-1]) / (2*deltay) - abs(v[0][0]) * (u[i][j+1] - 2*u[i][j] + u[i][j-1]) / (2*deltay)) - ¥
                ((1/rho) * (p[i+1][j]) / (deltax)) + ¥
                nu((u[i+1][j] - 2*u[i][j] + u[i][j-1]) / ((deltax)**2) + (u[i][j+1] - 2*u[i][j] + u[i][j-1]) / ((deltay)**2))
        print "T = " + str(T)
        print "u1 = " + str(u1)
        print "u2 = " + str(u2)
        print "u3 = " + str(u3)
        #print "u2 = " +str(u2)
        #print "u3 = " +str(u3)
        u += 1
        #u2 += 1
        #u3 += 1
    T += 1
print "End"
