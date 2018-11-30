# -*- coding: utf-8 -*-

import numpy as np
import math
import matplotlib.pyplot as plt
from numpy.random import *
from matplotlib import animation
import matplotlib.patches as patches
import os
import csv
import sys
value = sys.argv
os.chdir('/media/pascal/HD-GDU3/Tajima_backup/EHD/result')
os.mkdir(str(value[1]))
import time
start_time = time.time()

#物性値
print "物性値入力"
rho = input("rho(流体密度)[kg/m^3] = ")
nu = input("nu（動粘性係数）[m^2/s] = ")

#定数
print "定数入力"
H = input("H(流れ場y方向長さ)[m] = ")
L = input("L(流れ場x方向長さ)[m] = ")
T = input("T(移流時間)[s] = ")
deltaT = input("deltaT(時間刻み)[s] = ")
deltax = input("deltax(x方向要素間距離)[m] = ")
deltay = input("deltay(y方向要素間距離)[m] = ")
omega = input("omega(速度場補正緩和係数) = ")
gamma = input("gamma(圧力補正緩和係数) = ")
M = input("M(連続の式収束条件) = ")

#物性値、定数の出力
constant_list = [["rho",0],["nu",0],["H",0],["L",0],["T",0],["deltaT",0],["deltax",0],["deltay",0],["omega",0],["gamma", 0],["M",0]]
constant_list[0][1] = rho
constant_list[1][1] = nu
constant_list[2][1] = H
constant_list[3][1] = L
constant_list[4][1] = T
constant_list[5][1] = deltaT
constant_list[6][1] = deltax
constant_list[7][1] = deltay
constant_list[8][1] = omega
constant_list[9][1] = gamma
constant_list[10][1] = M
with open(os.path.join(str(value[1]),"constant.csv"), 'w') as file:
    writer = csv.writer(file, lineterminator = '\n')
    writer.writerows(constant_list)

#流れ場の座標(ms,n)
ms = int(L / deltax + 1)
n = int(H / deltay + 1)
print "ms :" + str(ms)
print "n :" + str(n)

#境界における初期値
u_BoundaryAD = 0.0#x方向速度[m/s]@inlet
v_BoundaryAD = 0.0#y方向速度[m/s]@inlet
p_BoundaryAD = 3#圧力[Pa]@inlet

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

#対流項CNVと粘性項DIF配列設定
CNVU = np.array([[0.0] * (n+1) for i in range(ms+1)])
CNVV = np.array([[0.0] * (n+1) for i in range(ms+1)])
DIFU = np.array([[0.0] * (n+1) for i in range(ms+1)])
DIFV = np.array([[0.0] * (n+1) for i in range(ms+1)])
#連続の式DIV・圧力補正量deltapの配列設定
DIV = np.array([[0.0] * (n+1) for i in range(ms+1)])

#圧力から速度場を決定する場合(全速度場指定)
def u_allsetting():
    i = 0
    j = 1
    while 0 <= i <= ms-1:
        while 1 <= j <= n-1:
            u_old[i][j] = 1.0 / (2*nu*rho) * (-1.0*(p_BoundaryBC-p_BoundaryAD)/ L) * (1.0*(j-0.5)*deltay) * (H-(1.0*(j-0.5)*deltay))
            j += 1
        j = 1
        i += 1
u_allsetting()

#全圧力場指定
def p_allsetting():
    i = 0
    j = 0
    while 0 <= i <= ms:
        while 0 <= j <= n:
            p_old[i][j] = p_BoundaryAD + ((p_BoundaryAD-p_BoundaryBC)/(ms-1)*(0.5-i))
            j += 1
        j = 0
        i += 1
p_allsetting()

#境界条件の設定
def boundary_condition():
    #BoundaryAD
    j = 0
    while 0 <= j <= n-1:
        u_old[0][j] = 1.2 / (2*nu*rho) * (-1.0*(p_BoundaryBC-p_BoundaryAD)/ L) * (1.0*(j-0.5)*deltay) * (H-(1.0*(j-0.5)*deltay))
        v_old[0][j] = 2.0 * v_BoundaryAD - v_old[1][j]
        j += 1
    #WallAB
    i = 0
    while 0 <= i <= ms-1:
        u_old[i][0] = 2.0 * u_WallAB - u_old[i][1]
        v_old[i][0] = v_WallAB
        i += 1
    #WallCD
    i = 0
    while 0 <= i <= ms-1:
        u_old[i][n] = 2.0 * u_WallCD - u_old[i][n-1]
        v_old[i][n-1] = v_WallCD
        i += 1
    #BoundaryBC
    j = 0
    while 0 <= j <= n-1:
        u_old[ms-1][j] = 1.2 / (2*nu*rho) * (-1.0*(p_BoundaryBC-p_BoundaryAD)/ L) * (1.0*(j-0.5)*deltay) * (H-(1.0*(j-0.5)*deltay))
        v_old[ms][j] = 2.0 * v_BoundaryBC - v_old[ms-1][j]
        j += 1
boundary_condition()

#初期におけるDIV
while 1 <= i <= ms-1:
    while 1 <= j <= n-1:
        DIV[i][j] = abs((u_old[i][j] - u_old[i-1][j])*1.0/deltax + (v_old[i][j] - v_old[i][j-1])*1.0/deltay)
        j += 1
    j = 1
    i += 1

t = 0
m = 1
#csvファイルで出力
def csvout():
    u_out = u_old.transpose()
    v_out = v_old.transpose()
    p_out = p_old.transpose()
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
    with open(os.path.join(str(value[1]),"DIV_(t="+str(t)+")"+".csv"), 'w') as file:
        writer = csv.writer(file, lineterminator = '\n')
        writer.writerows(DIV_out)
csvout()

#グラフを作成して保存する
X = np.array([[0.0] * (n+1) for i in range(ms+1)])
i = 0
j = 0
while 0 <= j <= n:
    while 0 <= i <= ms:
        X[i][j] = deltax * i
        i += 1
    i = 0
    j += 1
Y = np.array([[0.0] * (n+1) for i in range(ms+1)])
i = 0
j = 0
while 0 <= j <= n:
    while 0 <= i <= ms:
        Y[i][j] = deltay * j
        i += 1
    i = 0
    j += 1
X_out = X.transpose()
Y_out = Y.transpose()
def graph():
    #配列変換・設定
    u_out = u_old.transpose()
    v_out = v_old.transpose()
    p_out = p_old.transpose()
    velocity_out = np.sqrt(u_out**2+v_out**2)
    #ディレクトリ移動
    os.chdir(str(value[1]))
    #速度ベクトル作成
    plt.pcolor(X_out, Y_out, velocity_out)
    plt.colorbar()
    plt.quiver(X_out, Y_out, u_out, v_out, angles='xy', scale_units='xy', scale=np.max(velocity_out)*(1.0/(L*0.03/2)), headwidth=5, headlength=8, headaxislength=4)
    plt.axis('equal')
    if m % 10000 == 0:
        plt.title('velocity_vector(t='+str(t)+',m='+str(m)+')')
    else:
        plt.title('velocity_vector(t='+str(t)+')')
    plt.xlabel('x')
    plt.ylabel('y')
    plt.xlim(-1.0*L/10, 11.0*L/10)
    plt.ylim(-1.0*H/10, 11.0*H/10)
    plt.grid()
    plt.draw()
    ax = plt.axes()
    r = patches.Rectangle(xy=(0, -0.001), width=0.0005, height=0.0005, fc='y')
    ax.add_patch(r)
    if m % 10000 == 0:
        plt.savefig("velocity(t=" + str(t) +",m="+str(m)+ ").png", dpi=600)
    else:
        plt.savefig("velocity(t=" + str(t) + ").png", dpi=600)
    plt.cla()
    plt.clf()
    plt.close()
    #圧力分布作成
    plt.pcolor(X_out, Y_out, p_out)
    plt.colorbar()
    plt.axis('equal')
    if m % 10000 == 0:
        plt.title('pressure_distribution(t='+str(t)+',m='+str(m)+')')
    else:
        plt.title('pressure_distribution(t='+str(t)+')')
    plt.xlabel('x')
    plt.ylabel('y')
    plt.xlim(-1.0*L/10, 11.0*L/10)
    plt.ylim(-1.0*H/10, 11.0*H/10)
    plt.grid()
    if m % 10000 == 0:
        plt.savefig("pressure(t=" + str(t) +",m="+str(m)+ ").png", dpi=600)
    else:
        plt.savefig("pressure(t=" + str(t) + ").png", dpi=600)
    plt.cla()
    plt.clf()
    plt.close()
    os.chdir('../')
graph()

print "Calculation starts"
t = deltaT
while t <= T:
    print "t =" + str(t)
    #u_old,v_old仮値設定①粘性項・対流項配列の設定
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
    #u_old,v_old仮値設定②ナビエストークス方程式を解く
    i = 1
    j = 1
    while 1 <= i <= ms-2:
        while 1 <= j <= n-1:
            u_old[i][j] = u_old[i][j] + deltaT * (-(1.0/rho)*(p_old[i+1][j]-p_old[i][j])/deltax - CNVU[i][j] + DIFU[i][j])
            j += 1
        j = 1
        i += 1
    i = 1
    j = 1
    while 1 <= i <= ms-1 :
        while 1 <= j <= n-2:
            v_old[i][j] = v_old[i][j] + deltaT * (-(1.0/rho)*(p_old[i][j+1]-p_old[i][j])/deltay - CNVV[i][j] + DIFV[i][j])
            j += 1
        j = 1
        i += 1

    #圧力補正ループ
    m = 1
    Dmax = M + 1
    while Dmax > M :
        #print "配列DIV内の最大値DmaxがMより大きい場合ループに入る"
        print str(value[1])
        print "t = " + str(t)
        print "m = " + str(m)
        i = 1
        j = 1
        while 1 <= i <= ms-1:
            while 1 <= j <= n-1:
                deltap = - rho * 1.0 / (2*deltaT) * (deltax*deltay) / (deltax**2+deltay**2) * (deltay * (u_old[i][j]-u_old[i-1][j]) + deltax * (v_old[i][j]-v_old[i][j-1]))
                p_old[i][j] = p_old[i][j] + 1.0 * gamma * deltap
                if i <= ms-2:
                    u_old[i][j] = u_old[i][j] + 1.0 * omega * (1.0/rho) * (deltaT/deltax) * deltap
                if j <= n-2:
                    v_old[i][j] = v_old[i][j] + 1.0 * omega * (1.0/rho) * (deltaT/deltay) * deltap
                if 2 <= i:
                    u_old[i-1][j] = u_old[i-1][j] - 1.0 * omega * (1.0/rho) * (deltaT/deltax) * deltap
                if 2 <= j:
                    v_old[i][j-1] = v_old[i][j-1] - 1.0 * omega * (1.0/rho) * (deltaT/deltay) * deltap
                j += 1
            j = 1
            i += 1

        boundary_condition()

        #連続の式の収束条件
        i = 1
        j = 1
        while 1 <= i <= ms-1:
            while 1 <= j <= n-1:
                DIV[i][j] = abs((u_old[i][j] - u_old[i-1][j])*1.0/deltax + (v_old[i][j] - v_old[i][j-1])*1.0/deltay)
                j += 1
            j = 1
            i += 1
        #---通常モード---
        Dmax = np.max(DIV)
        #---↑---
        #---↓Dmaxの発散を検知したら補正ループを終了する↓---
        #Dmax_new = np.max(DIV)
        #if Dmax_new > Dmax:#発散：自動的に計算終了
            #print  "Divergence in Dmax = " + str(Dmax_new)
            #Dmax = 0
        #else:#収束：Dmaxを更新して計算続行
            #if m > 500:
                #print  "stop at m = " + str(m) +"(Dmax = " + str(Dmax_new) + ")"
                #Dmax = 0
            #else:
                #Dmax = Dmax_new
                #print "Converging now"
        #---↑---
        #---↓圧力補正計算を回数で制御する↓---
        #Dmax_new = np.max(DIV)
        #if m >200:
            #print  "stop at m = " + str(m) +"(Dmax = " + str(Dmax_new) + ")"
            #Dmax = 0
        #else:
            #Dmax = Dmax_new
            #print "Converging now"
        #---↑---
        print "Dmax = " + str(Dmax)
        print "deltap = " + str(deltap)
        process_time = time.time() - start_time
        print "process_time = " + str(process_time)
        if m % 10000 == 0:
            csvout()
            graph()
        m += 1
        #Dmax = 0#強制的ループ終了用
    #csvファイルとグラフを出力
    if t == deltaT or int(t/deltaT) % 20 == 0:
        csvout()
        graph()

    #時間を進める
    t = t + deltaT
print "Calculation ends"
