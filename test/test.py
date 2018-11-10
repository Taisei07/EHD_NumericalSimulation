# -*- coding: utf-8 -*-

import numpy as np
import matplotlib.pyplot as plt
from numpy.random import *
import os

#圧力補正量を設定する

#物性値
rho = 1000.0 #密度[kg/m^3]
nu = 0.000001#動粘度[m^2/s]
#定数
ms = 800#x方向の要素数
n = 40#y方向の要素数
T = 5#移流時間
deltaT = 1.0#時間刻み
deltax = 0.0025#x方向の要素間距離
deltay = 0.0025#y方向の要素間距離
omega = 0.5#圧力補正における緩和係数
M = 0.01#連続の式の収束条件
#初期条件
u_boundaryAH = 0.9#入り口x方向速度[m/s]
u_boundaryFG = 0.0#出口x方向速度[m/s]
u_wallAF = 0.0#wallAFx方向速度[m/s]
u_wallGH = 0.0#wallGHx方向速度[m/s]
p_boundaryAH = 0#入り口圧力[Pa]
p_boundaryFG = 0.0#出口圧力[Pa]
p_wallAF = 0.0
p_wallGH = 0.0

#格子定義点(oldはtimestep=n, newはtimestep=n+1)
u_old = np.array([[0.0] * n for i in range(ms)])#ms:x方向,n:y方向
v_old = np.array([[0.0] * n for i in range(ms)])
p_old = np.array([[0.0] * n for i in range(ms)])
u_new = np.array([[0.0] * n for i in range(ms)])
v_new = np.array([[0.0] * n for i in range(ms)])
p_new = np.array([[0.0] * n for i in range(ms)])
#u_old,v_oldには初期値（乱数）を与える必要がある
#i = 0
#j = 0
#while 0 <= i <= ms-1:
    #while 0 <= j <= n-1:
        #u_old[i][j] = 1.0 * randint(1,10)
        #j += 1
    #j = 0
    #i += 1
#while 0 <= i <= ms-1:
    #while 0 <= j <= n-1:
        #v_old[i][j] = 1.0 * randint(1,10)
        #j += 1
    #j = 0
    #i += 1
#対流項CNVと粘性項DIF配列設定
CNVU = np.array([[0.0] * n for i in range(ms)])
CNVV = np.array([[0.0] * n for i in range(ms)])
DIFU = np.array([[0.0] * n for i in range(ms)])
DIFV = np.array([[0.0] * n for i in range(ms)])
#連続の式用のDIV配列設定
DIV = np.array([[0.0] * n for i in range(ms)])
#圧力補正量deltap配列設定
deltap = np.array([[0.0] * n for i in range(ms)])
#初期条件の設定
#boundaryAH(inlet)境界
i = 0
j = 1
while j <= n-2:
    u_old[0][j] = u_boundaryAH
    p_old[0][j] = 0
    p_old[1][j] = 2*p_boundaryAH-p_old[0][j]
    j += 1
print "inlet初期条件完了"
#boundaryFG(outlet)境界
i = 0
j = 1
while j <= n-2:
    u_old[ms-1][j] = u_boundaryFG
    p_old[ms-1][j] = 0
    p_old[ms-2][j] = 2*p_boundaryFG-p_old[ms-1][j]
    j += 1
print "outlet初期条件完了"
#wallAF境界
i = 2
j = 0
while i <= ms-2:
    u_old[i][0] = 0
    u_old[i][1] = 2*u_wallAF-u_old[i][0]
    p_old[i][0] = 0
    p_old[i][1] = 2*p_wallAF-p_old[i][0]
    i += 1
print "wall_AF初期条件完了"
#wallGH境界
i = 2
j = 0
while i <= ms-2:
    u_old[i][n-1] = 0
    u_old[i][n-2] = 2*u_wallGH-u_old[i][n-1]
    p_old[i][n-1] = 0
    p_old[i][n-2] = 2*p_wallGH-p_old[i][n-1]
    i += 1
print "wall_GH初期条件完了"

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
t = 0
while t <= T:
    print "t =" + str(t)
    #u_old,v_oldの仮値を設定①粘性項・対流項配列の設定
    i = 1
    j = 1
    while 0 < i < ms-1:
        while 0 < j < n-1:
            CNVU[i][j] = (u_old[i][j] * (u_old[i+1][j]-u_old[i-1][j]) * 1.0 /(2*deltax) - abs(u_old[i][j]) * (u_old[i+1][j]-2*u_old[i][j]+u_old[i-1][j]) * 1.0 / (2*deltax)) - (v_old[i][j]*(u_old[i][j+1] - u_old[i][j-1]) * 1.0 / (2 * deltay) - abs(v_old[i][j])*(u_old[i][j+1]-2*u_old[i][j]+u_old[i][j-1]) * 1.0 / deltay)
            CNVV[i][j] = (u_old[i][j] * (v_old[i+1][j]-v_old[i-1][j]) * 1.0 /(2*deltax) - abs(u_old[i][j]) * (v_old[i+1][j]-2*v_old[i][j]+v_old[i-1][j]) * 1.0 / (2*deltax)) - (v_old[i][j]*(v_old[i][j+1] - v_old[i][j-1]) * 1.0 / (2 * deltay) - abs(v_old[i][j])*(v_old[i][j+1]-2*v_old[i][j]+v_old[i][j-1]) * 1.0 / deltay)
            DIFU[i][j] = nu * ((u_old[i+1][j]-2*u_old[i][j]+u_old[i-1][j]) * 1.0 / ((deltax)**2) + (u_old[i][j+1]-2*u_old[i][j]+u_old[i][j-1]) * 1.0 / ((deltay)**2))
            DIFV[i][j] = nu * ((v_old[i+1][j]-2*v_old[i][j]+v_old[i-1][j]) * 1.0 / ((deltax)**2) + (v_old[i][j+1]-2*v_old[i][j]+v_old[i][j-1]) * 1.0 / ((deltay)**2))
            j += 1
        j = 1
        i += 1
    #u_old,v_oldの仮値を設定②ナビエストークス方程式を解く
    i = 1
    j = 1
    while 0 < i < ms-1:
        while 0 < j < n-1:
            u_old[i][j] = deltaT * (u_old[i][j] - (1.0/rho)*(p_old[i+1][j]-p_old[i][j])/deltax - CNVU[i][j] + DIFU[i][j])
            v_old[i][j] = deltaT * (v_old[i][j] - (1.0/rho)*(p_old[i][j+1]-p_old[i][j])/deltay - CNVV[i][j] + DIFV[i][j])
            j += 1
        j = 1
        i += 1
        #print "x方向にループ計算中"
    print u_old
    print v_old
    print "u,vの仮値設定完了"
    #u,vの仮値におけるDIVを設定
    i = 1
    j = 1
    while 0 < i < ms-1:
        while 0 < j < n-1:
            DIV[i][j] =  abs((u_old[i][j] - u_old[i-1][j])*1.0/deltax + (v_old[i][j] - v_old[i][j-1])*1.0/deltay)
            j += 1
        j = 1
        i += 1
    print "仮値のDIV"
    print DIV
    Dmax = np.max(DIV)
    print "u,vの仮値におけるDIVの配列設定完了"

    #圧力補正量の導出
    i = 1
    j = 1
    m = 1#反復回数
    while Dmax > M :
        print "配列DIV内の最大値DmaxがMより大きい場合ループに入る"
        print "Dmax = " + str(Dmax)
        while 0 < i < ms-1:
            #print "全格子点に対する圧力補正量deltapを設定"
            while 0 < j < n-1:
                deltap[i][j] = -rho*1.0/(2*deltaT)*(deltax * deltay)/(deltax**2 + deltay**2)*(deltay*(u_old[i][j]-u_old[i-1][j])+deltax*(v_old[i][j]-v_old[i][j-1]))
                j += 1
            j = 1
            i += 1
        print "m=" + str(m) + "の時の圧力補正量deltapの設定完了"
        i = 1
        j = 1
        while 0 < i < ms-1:
            #print "設定したdeltapを適用して、速度場を更新"
            while 0 < j < n-1:
                p_old[i][j] = p_old[i][j] + (omega * deltap[i][j])
                u_old[i][j] = u_old[i][j] + 1.0 * omega * (1/rho) * (deltaT/deltax) * deltap[i][j]
                v_old[i][j] = v_old[i][j] + 1.0 * omega * (1/rho) * (deltaT/deltay) * deltap[i][j]
                j += 1
            j = 1
            i += 1
        i = 1
        j = 1
        while 0 < i < ms-1:
            while 0 < j < n-1:
                u_old[i-1][j] = u_old[i-1][j] - 1.0 * omega * (1/rho) * (deltaT/deltax) * deltap[i][j]
                v_old[i][j-1] = u_old[i][j-1] - 1.0 * omega * (1/rho) * (deltaT/deltay) * deltap[i][j]
                j += 1
            j = 1
            i += 1
        #境界条件
        #boundaryAH(inlet)境界
        i = 0
        j = 0
        while j <= n-1:
            u_old[0][j] = u_boundaryAH
            p_old[0][j] = 0
            p_old[1][j] = 2*p_boundaryAH-p_old[0][j]
            j += 1
        #boundaryFG(outlet)境界
        i = 0
        j = 0
        while j <= n-1:
            u_old[ms-1][j] = u_boundaryFG
            p_old[ms-1][j] = 0
            p_old[ms-2][j] = 2*p_boundaryFG-p_old[ms-1][j]
            j += 1
        #wallAF境界
        i = 0
        j = 0
        while i <= ms-1:
            u_old[i][0] = 0
            u_old[i][1] = 2*u_wallAF-u_old[i][0]
            p_old[i][0] = 0
            p_old[i][1] = 2*p_wallAF-p_old[i][0]
            i += 1
        #wallGH境界
        i = 0
        j = 0
        while i <= ms-1:
            u_old[i][n-1] = 0
            u_old[i][n-2] = 2*u_wallGH-u_old[i][n-1]
            p_old[i][n-1] = 0
            p_old[i][n-2] = 2*p_wallGH-p_old[i][n-1]
            i += 1
        print "m=" + str(m) + "の時の速度場u,vの設定完了"
        print "p_old"
        print p_old
        print "u_old"
        print u_old
        print "v_old"
        print v_old
        i = 1
        j = 1
        while 0 < i < ms-1:
            while 0 < j < n-1:
                DIV[i][j] =  abs((u_old[i][j] - u_old[i-1][j])*1.0/deltax + (v_old[i][j] - v_old[i][j-1])*1.0/deltay)
                j += 1
                print DIV[i][j]
            j = 1
            i += 1
        i = 1
        j = 1
        print "DIV"
        print DIV
        Dmax = np.max(DIV)
        m += 1
        #Dmax = 0#強制的ループ終了用
    print "全格子点が連続の式の収束条件を満たす"

    print "u_old[50][18] = " + str(u_old[50][18]) + "u_old[19][50] = " + str(u_old[51][18])
    print "v_old[50][18] = " + str(v_old[50][18]) + "v_old[19][50] = " + str(v_old[51][18])
    print "v_old[40][10] = " + str(v_old[40][10]) + "v_old[10][41] = " + str(v_old[41][10])
    print p_old

    #csvファイルで出力
    u_out = u_old.transpose()
    v_out = v_old.transpose()
    p_out = p_old.transpose()
    deltap_out = deltap.transpose()
    DIV_out = DIV.transpose()
    import csv
    with open("u_(t="+str(t)+")"+".csv", 'w') as file:
        writer = csv.writer(file, lineterminator = '\n')
        writer.writerows(u_out)
    with open("v_(t="+str(t)+")"+".csv", 'w') as file:
        writer = csv.writer(file, lineterminator = '\n')
        writer.writerows(v_out)
    with open("p_(t="+str(t)+")"+".csv", 'w') as file:
        writer = csv.writer(file, lineterminator = '\n')
        writer.writerows(p_out)
    with open("deltap_(t="+str(t)+")"+".csv", 'w') as file:
        writer = csv.writer(file, lineterminator = '\n')
        writer.writerows(deltap_out)
    with open("DIV_(t="+str(t)+")"+".csv", 'w') as file:
        writer = csv.writer(file, lineterminator = '\n')
        writer.writerows(DIV_out)

    #時間を進める
    t = t + deltaT
print "計算終了"
