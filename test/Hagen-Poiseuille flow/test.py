# -*- coding: utf-8 -*-

import numpy as np
import matplotlib.pyplot as plt
from numpy.random import *
import os
import sys
value = sys.argv
os.mkdir(str(value[1]))

#圧力補正量を設定する

#物性値
rho = 1000.0 #密度[kg/m^3]
nu = 0.000001#動粘度[m^2/s]
#定数
H =0.1#流れ場のy方向の長さ
L =2.0#流れ場のx方向の長さ
T = 0.00001#移流時間
deltaT = 0.000001#時間刻み
deltax = 0.0025#x方向の要素間距離
deltay = 0.0025#y方向の要素間距離
omega = 0.5#圧力補正における緩和係数
M = 0.00001#連続の式の収束条件

#流れ場の座標(ms,n)
ms = int(L / deltax + 1)
n = int(H / deltay + 1)
print "ms :" + str(ms)
print "n :" + str(n)

#初期条件
u_BoundaryAD = 0.0#x方向速度[m/s]@inlet
v_BoundaryAD = 0.0#y方向速度[m/s]@inlet
p_BoundaryAD = 100.0#圧力{Pa}@inlet

u_WallAB = 0.0#x方向速度[m/s]@wallAB
v_WallAB = 0.0#x方向速度[m/s]@wallAB
p_WallAB = 0.0#圧力[Pa]@wallAB

u_BoundaryBC = 0.0#x方向速度[m/s]@outlet
v_BoundaryBC = 0.0#y方向速度[m/s]@outlet
p_BoundaryBC = 0.0#圧力[Pa]@outlet

u_WallCD = 0.0#x方向速度[m/s]@wallCD
v_WallCD = 0.0#y方向速度[m/s]@wallCD
p_WallCD = 0.0#圧力[Pa]@wallCD

#格子定義点(oldはtimestep=n, newはtimestep=n+1)　＊座標(ms,n)に1を足している
u_old = np.array([[0.0] * (n+1) for i in range(ms+1)])#ms:x方向,n:y方向
v_old = np.array([[0.0] * (n+1) for i in range(ms+1)])
p_old = np.array([[0.0] * (n+1) for i in range(ms+1)])
u_new = np.array([[0.0] * (n+1) for i in range(ms+1)])
v_new = np.array([[0.0] * (n+1) for i in range(ms+1)])
p_new = np.array([[0.0] * (n+1) for i in range(ms+1)])

#対流項CNVと粘性項DIF配列設定
CNVU = np.array([[0.0] * (n+1) for i in range(ms+1)])
CNVV = np.array([[0.0] * (n+1) for i in range(ms+1)])
DIFU = np.array([[0.0] * (n+1) for i in range(ms+1)])
DIFV = np.array([[0.0] * (n+1) for i in range(ms+1)])
#連続の式DIV・圧力補正量deltapの配列設定
DIV = np.array([[0.0] * (n+1) for i in range(ms+1)])

#圧力から速度場を決定する場合
j = 1
while 1 <= j <= n-1:
    u_old[0][j] = 1.0 / (2*nu*rho) * (-1.0*(p_BoundaryBC-p_BoundaryAD)/ L) * (1.0*(j-0.5)*deltay) * (H-(1.0*(j-0.5)*deltay))
    j += 1
print np.max(u_old)

#流れ場及び圧力場に放物型の解を与える(今回だけ)
i = 1
j = 1
while 1 <= i <= ms-1:
    while 1 <= j <= n-1:
        u_old[i][j] = 1.0 / (2*nu*rho) * (-1.0*(p_BoundaryBC-p_BoundaryAD)/ L) * (1.0*(j-0.5)*deltay) * (H-(1.0*(j-0.5)*deltay))
        j += 1
    j = 1
    i += 1
print u_old
i = 0
j = 0
while 0 <= i <= ms:
    while 0 <= j <= n:
        p_old[i][j] = 100.0625 - i * 0.125
        j += 1
    j = 0
    i += 1
print p_old

#境界条件の設定
#BoundaryAD
j = 0
while 0 <= j <= n-1:
    #u_old[0][j] = u_BoundaryAD
    v_old[0][j] = 2 * v_BoundaryAD - v_old[1][j]
    j += 1
#WallAB
i = 0
while 0 <= i <= ms-1:
    u_old[i][0] = 2 * u_WallAB - u_old[i][1]
    v_old[i][0] = v_WallAB
    i += 1
#WallCD
i = 0
while 0 <= i <= ms-1:
    u_old[i][n] = 2 * u_WallCD - u_old[i][n-1]
    v_old[i][n-1] = v_WallCD
    i += 1
#BoundaryBC
j = 0
while 0 <= j <= n-1:
    #u_old[ms-1][j] = u_BoundaryBC
    v_old[ms][j] = 2 * v_BoundaryBC - v_old[ms-1][j]
    j += 1


#初期条件の確認
print u_old
print v_old
print p_old

print "ここから計算開始"
t = 0
while t <= T:
    print "t =" + str(t)
    #u_old,v_oldの仮値を設定①粘性項・対流項配列の設定
    i = 1
    j = 1
    while 1 <= i <= ms-1:
        while 1 <= j <= n-1:
            CNVU[i][j] = (u_old[i][j] * (u_old[i+1][j]-u_old[i-1][j]) * 1.0 /(2*deltax) - abs(u_old[i][j]) * (u_old[i+1][j]-2*u_old[i][j]+u_old[i-1][j]) * 1.0 / (2*deltax)) - (v_old[i][j]*(u_old[i][j+1] - u_old[i][j-1]) * 1.0 / (2 * deltay) - abs(v_old[i][j])*(u_old[i][j+1]-2*u_old[i][j]+u_old[i][j-1]) * 1.0 / (2 * deltay))
            CNVV[i][j] = (u_old[i][j] * (v_old[i+1][j]-v_old[i-1][j]) * 1.0 /(2*deltax) - abs(u_old[i][j]) * (v_old[i+1][j]-2*v_old[i][j]+v_old[i-1][j]) * 1.0 / (2*deltax)) - (v_old[i][j]*(v_old[i][j+1] - v_old[i][j-1]) * 1.0 / (2 * deltay) - abs(v_old[i][j])*(v_old[i][j+1]-2*v_old[i][j]+v_old[i][j-1]) * 1.0 / (2 * deltay))
            DIFU[i][j] = nu * ((u_old[i+1][j]-2*u_old[i][j]+u_old[i-1][j]) * 1.0 / ((deltax)**2) + (u_old[i][j+1]-2*u_old[i][j]+u_old[i][j-1]) * 1.0 / ((deltay)**2))
            DIFV[i][j] = nu * ((v_old[i+1][j]-2*v_old[i][j]+v_old[i-1][j]) * 1.0 / ((deltax)**2) + (v_old[i][j+1]-2*v_old[i][j]+v_old[i][j-1]) * 1.0 / ((deltay)**2))
            j += 1
        j = 1
        i += 1
    #u_old,v_oldの仮値を設定②ナビエストークス方程式を解く
    i = 1
    j = 1
    while 1 <= i <= ms-2:
        while 1 <= j <= n-1:
            u_old[i][j] = u_old[i][j] + deltaT * (-(1.0/rho)*(p_old[i+1][j]-p_old[i][j])/deltax - CNVU[i][j] + DIFU[i][j])
            j += 1
        j = 1
        i += 1
    while 1 <= i <= ms-1 :
        while 1 <= j <= n-2:
            v_old[i][j] = v_old[i][j] + deltaT * (-(1.0/rho)*(p_old[i][j+1]-p_old[i][j])/deltay - CNVV[i][j] + DIFV[i][j])
            j += 1
        j = 1
        i += 1
    print u_old
    print v_old
    print "u,vの仮値設定完了"
    #print str(u_old[1][1])
    #print str(u_old[0][1])

    #u,vの仮値におけるDIVを設定
    i = 1
    j = 1
    while 1 <= i <= ms-1:
        while 1 <= j <= n-1:
            DIV[i][j] = abs((u_old[i][j] - u_old[i-1][j])*1.0/deltax) + abs((v_old[i][j] - v_old[i][j-1])*1.0/deltay)
            j += 1
        j = 1
        i += 1
    print "仮値のDIV"
    print DIV
    Dmax = np.max(DIV)
    print Dmax
    print "u,vの仮値におけるDIVの配列設定完了"

    #圧力補正量の導出
    m = 1#反復回数
    while Dmax > M :
        #print "配列DIV内の最大値DmaxがMより大きい場合ループに入る"
        print "m=" + str(m)
        print "Dmax = " + str(Dmax)
        i = 1
        j = 1
        while 1 <= i <= ms-1:
            while 1 <= j <= n-1:
                deltap = -rho*1.0/(2*deltaT)*(deltax * deltay)/(deltax**2 + deltay**2)*(deltay*(u_old[i][j]-u_old[i-1][j])+deltax*(v_old[i][j]-v_old[i][j-1]))
                p_old[i][j] = p_old[i][j] + (omega * deltap)
                if i <= ms-2:
                    u_old[i][j] = u_old[i][j] + 1.0 * omega * (1/rho) * (deltaT/deltax) * deltap
                if j <= n-2:
                    v_old[i][j] = v_old[i][j] + 1.0 * omega * (1/rho) * (deltaT/deltay) * deltap
                if 2 <= i:
                    u_old[i-1][j] = u_old[i-1][j] - 1.0 * omega * (1/rho) * (deltaT/deltax) * deltap
                if 2 <= j:
                    v_old[i][j-1] = v_old[i][j-1] - 1.0 * omega * (1/rho) * (deltaT/deltay) * deltap
                j += 1
            j = 1
            i += 1
        j = 1
        print u_old
        print v_old
        print p_old
        #境界条件の適用
        #BoundaryAD
        j = 1
        while 1 <= j <= n-2:
            v_old[0][j] = 2 * v_BoundaryAD - v_old[1][j]
            j += 1
        #WallAB
        i = 1
        while 1 <= i <= ms-2:
            u_old[i][0] = 2 * u_WallAB - u_old[i][1]
            i += 1
        #WallCD
        i = 1
        while 1 <= i <= ms-2:
            u_old[i][n] = 2 * u_WallCD - u_old[i][n-1]
            i += 1
        #BoundaryBC
        j = 1
        while 1 <= j <= n-2:
            v_old[ms][j] = 2 * v_BoundaryBC - v_old[ms-1][j]
            j += 1
        #連続の式の収束条件
        i = 1
        j = 1
        while 1 <= i <= ms-1:
            while 1 <= j <= n-1:
                DIV[i][j] = abs((u_old[i][j] - u_old[i-1][j])*1.0/deltax) + abs((v_old[i][j] - v_old[i][j-1])*1.0/deltay)
                j += 1
            j = 1
            i += 1
        print DIV
        #---通常モード---
        #Dmax = np.max(DIV)
        #---↑---
        #---↓Dmaxの発散を検知したら補正ループを終了するプログラミング↓---
        Dmax_new = np.max(DIV)
        if Dmax_new > Dmax:#発散：自動的に計算終了
            print  "Divergence in Dmax = " + str(Dmax_new)
            Dmax = 0
        else:#収束：Dmaxを更新して計算続行
            if m > 200:
                print  "stop at m = " + str(m) +"(Dmax = " + str(Dmax_new) + ")"
                Dmax = 0
            else:
                Dmax = Dmax_new
                print "Converging now"
        #---↑---
        m += 1
        #Dmax = 0#強制的ループ終了用
    #csvファイルで出力
    u_out = u_old.transpose()
    v_out = v_old.transpose()
    p_out = p_old.transpose()
    #deltap_out = deltap.transpose()
    DIV_out = DIV.transpose()
    import csv
    with open(os.path.join(str(value[1]),"u_(t="+str(t)+")"+".csv"), 'w') as file:
        writer = csv.writer(file, lineterminator = '\n')
        writer.writerows(u_out)
    with open(os.path.join(str(value[1]),"v_(t="+str(t)+")"+".csv"), 'w') as file:
        writer = csv.writer(file, lineterminator = '\n')
        writer.writerows(v_out)
    with open(os.path.join(str(value[1]),"p_(t="+str(t)+")"+".csv"), 'w') as file:
        writer = csv.writer(file, lineterminator = '\n')
        writer.writerows(p_out)
    #with open("deltap_(t="+str(t)+")"+".csv", 'w') as file:
        #writer = csv.writer(file, lineterminator = '\n')
        #writer.writerows(deltap_out)
    with open(os.path.join(str(value[1]),"DIV_(t="+str(t)+")"+".csv"), 'w') as file:
        writer = csv.writer(file, lineterminator = '\n')
        writer.writerows(DIV_out)

    #時間を進める
    t = t + deltaT
print str(value[1]) + "計算終了"
