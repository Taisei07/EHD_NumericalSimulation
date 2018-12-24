# -*- coding: utf-8 -*-

import numpy as np
import math
import matplotlib.pyplot as plt
import matplotlib.patches as patches
import json
import requests
import os
import sys
sys.path.append('../../../')
from module.LineAPI import line_notify_token
import csv
value = sys.argv
os.chdir('/media/pascal/HD-GDU3/Tajima_backup/EHD/result')
os.mkdir(str(value[1]))
import time
start_time = time.time()
import requests

#Line通知
def LineMessage():
    process_time = time.time() - start_time
    line_notify_api = 'https://notify-api.line.me/api/notify'
    headers = {'Authorization': 'Bearer ' + line_notify_token}
    # メッセージ
    payload = {'message': str(value[1]) + '\n' + 'process_time : ' + str(process_time)}
    requests.post(line_notify_api, data=payload, headers=headers)

def LineFigure(A):
    os.chdir(str(value[1]))
    line_notify_api = 'https://notify-api.line.me/api/notify'
    headers = {'Authorization': 'Bearer ' + line_notify_token}
    # メッセージ
    payload = {'message': str(A)}
    files = {"imageFile": open(A, 'rb')}
    requests.post(line_notify_api, data=payload, headers=headers, files=files)
    os.chdir('../')

#物性値
print "物性値入力"
rho = 1000#input("rho(流体密度)[kg/m^3] = ")
nu = 0.000001#input("nu（動粘性係数）[m^2/s] = ")

#定数
print "定数入力"
H = 0.002#input("H(流れ場y方向長さ)[m] = ")
L = 0.006#input("L(流れ場x方向長さ)[m] = ")
T = 0.5#input("T(移流時間)[s] = ")
deltaT = 0.00001#input("deltaT(時間刻み)[s] = ")
deltax = 0.0001#input("deltax(x方向要素間距離)[m] = ")
deltay = 0.0001#input("deltay(y方向要素間距離)[m] = ")
omega = 0.5#input("omega(緩和係数) = ")
M1 = 0.000001#input("M1(連続の式収束条件) = ")
M2 = 0.00001

#物性値、定数の出力
constant_list = [["rho",rho],["nu",nu],["H",H],["L",L],["T",T],["deltaT",deltaT],["deltax",deltax],["deltay",deltay],["omega",omega],["M1",M1]]
with open(os.path.join(str(value[1]),"constant.csv"), 'w') as file:
    writer = csv.writer(file, lineterminator = '\n')
    writer.writerows(constant_list)

#流れ場の座標(ms,n)
ms = int(L / deltax + 1)
n = int(H / deltay + 1)
print "ms :" + str(ms)
print "n :" + str(n)

#初期条件
u_BoundaryAD = 0.0#x方向速度[m/s]@inlet
v_BoundaryAD = 0.0#y方向速度[m/s]@inlet
p_BoundaryAD = 50.0#圧力[Pa]@inlet

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
deltap = np.array([[0.0] * (n+1) for i in range(ms+1)])

#圧力から速度場を決定する場合
if p_BoundaryAD != p_BoundaryBC:
    #速度場設定
    u_old[0:ms, 1:n] = [1.0 / (2*nu*rho) * (-1.0*(p_BoundaryBC-p_BoundaryAD)/ L) * (1.0*(i-0.5)*deltay) * (H-(1.0*(i-0.5)*deltay)) for i in range(1,n)]
    #圧力場指定
    p_old[0:ms+1, 0:n+1] = np.array([[p_BoundaryAD + ((p_BoundaryAD-p_BoundaryBC)/(ms-1)*(0.5-i))] for i in range(ms+1)])
#境界条件の設定
def boundary_condition():
    #BoundaryAD
    u_old[0:1, 0:n] = u_old[1:2, 0:n]
    v_old[0:1, 0:n] = v_old[1:2, 0:n]
    p_old[0:1, 0:n] = p_BoundaryAD
    #WallAB
    u_old[0:ms, 0:1] = -u_old[0:ms, 1:2]
    v_old[0:ms, 0:1] = v_WallAB
    p_old[0:ms, 0:1] = p_old[0:ms, 1:2]
    #WallCD
    u_old[0:ms, n:n+1] = -u_old[0:ms, n-1:n]
    v_old[0:ms, n-1:n] = v_WallCD
    p_old[0:ms, n:n+1] = p_old[0:ms, n-1:n]
    #BoundaryBC
    u_old[ms-1:ms, 0:n] = u_old[ms-2:ms-1, 0:n]
    v_old[ms:ms+1, 0:n] = v_old[ms-1:ms, 0:n]
    p_old[ms:ms+1, 0:n] = p_BoundaryBC
boundary_condition()

#初期におけるDIV
DIV[1:ms,1:n] = np.array([[abs((u_old[i][j] - u_old[i-1][j])*1.0/deltax + (v_old[i][j] - v_old[i][j-1])*1.0/deltay) for j in range(1,n)] for i in range(1,ms)])
#初期値を出力
t = 0
m1 = 1
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
X[0:ms+1, 0:n+1] = np.array([[deltax*(i-0.5)] for i in range(ms+1)])
Y = np.array([[0.0] * (n+1) for i in range(ms+1)])
Y[0:ms+1, 0:n+1] = np.array([deltay*(i-0.5) for i in range(n+1)])
X_out = X.transpose()
Y_out = Y.transpose()

def graph():
    #配列変換・設定
    u_out = u_old.transpose()
    v_out = v_old.transpose()
    p_out = p_old.transpose()
    velocity_out = np.sqrt(u_out**2+v_out**2)
    u_out[0:1, :] = 0
    u_out[n:n+1, :] = 0
    v_out[0:1, :] = 0
    v_out[n:n+1, :] = 0
    u_out[:, 0:1] = 0
    u_out[:, ms:ms+1] = 0
    v_out[:, 0:1] = 0
    v_out[:, ms:ms+1] = 0
    #ディレクトリ移動
    os.chdir(str(value[1]))
    #速度ベクトル作成
    plt.quiver(X_out, Y_out, u_out, v_out, angles='xy', scale_units='xy', scale=np.max(velocity_out)*(0.5/(L*0.03/2)), width=0.002, headwidth=7, headlength=10, headaxislength=5)
    plt.quiver(5.0/6*L, -1.0/5*H, np.max(velocity_out)*3, 0, angles='xy', scale_units='xy', scale=np.max(velocity_out)*(0.5/(L*0.03/2)), width=0.002, headwidth=7, headlength=10, headaxislength=5)
    plt.text(5.0/6*L, -2.0/5*H, str('{:.2E}'.format(np.max(velocity_out))) + '[m/s]')
    plt.axis('equal')
    if m1 % 10000 == 0:
        plt.title('velocity_vector(t='+str(t)+',m1='+str(m1)+')')
    else:
        plt.title('velocity_vector(t='+str(t)+')')
    plt.xlabel('x')
    plt.ylabel('y')
    plt.xlim(-1.0*L/10, 11.0*L/10)
    plt.ylim(-1.0*H/10, 11.0*H/10)
    plt.grid()
    ax = plt.axes()
    r = patches.Rectangle(xy=(0, 0), width=L, height=H, ec='#000000', fill=False)
    ax.add_patch(r)
    plt.draw()
    if m1 % 10000 == 0:
        plt.savefig("velocity(t=" + str(t) +",m1="+str(m1)+ ").png", dpi=600)
    else:
        plt.savefig("velocity(t=" + str(t) + ").png", dpi=600)
    plt.cla()
    plt.clf()
    plt.close()
    #圧力分布作成
    plt.pcolor(X_out, Y_out, p_out)
    plt.colorbar()
    plt.axis('equal')
    if m1 % 10000 == 0:
        plt.title('pressure_distribution(t='+str(t)+',m1='+str(m1)+')')
    else:
        plt.title('pressure_distribution(t='+str(t)+')')
    plt.xlabel('x')
    plt.ylabel('y')
    plt.xlim(-1.0*L/10, 11.0*L/10)
    plt.ylim(-1.0*H/10, 11.0*H/10)
    plt.grid()
    ax = plt.axes()
    r = patches.Rectangle(xy=(0, 0), width=L, height=H, ec='#000000', fill=False)
    ax.add_patch(r)
    if m1 % 10000 == 0:
        plt.savefig("pressure(t=" + str(t) +",m1="+str(m1)+ ").png", dpi=600)
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
    CNVU[1:ms, 1:n] = np.array([[(u_old[i][j] * (u_old[i+1][j]-u_old[i-1][j]) * 1.0 /(2*deltax) - abs(u_old[i][j]) * (u_old[i+1][j]-2*u_old[i][j]+u_old[i-1][j]) * 1.0 / (2*deltax)) - (v_old[i][j]*(u_old[i][j+1] - u_old[i][j-1]) * 1.0 / (2 * deltay) - abs(v_old[i][j])*(u_old[i][j+1]-2*u_old[i][j]+u_old[i][j-1]) * 1.0 / (2 * deltay)) for j in range(1, n)] for i in range(1,ms)])
    CNVV[1:ms, 1:n] = np.array([[(u_old[i][j] * (v_old[i+1][j]-v_old[i-1][j]) * 1.0 /(2*deltax) - abs(u_old[i][j]) * (v_old[i+1][j]-2*v_old[i][j]+v_old[i-1][j]) * 1.0 / (2*deltax)) - (v_old[i][j]*(v_old[i][j+1] - v_old[i][j-1]) * 1.0 / (2 * deltay) - abs(v_old[i][j])*(v_old[i][j+1]-2*v_old[i][j]+v_old[i][j-1]) * 1.0 / (2 * deltay)) for j in range(1, n)] for i in range(1,ms)])
    DIFU[1:ms, 1:n] = np.array([[nu * ((u_old[i+1][j]-2*u_old[i][j]+u_old[i-1][j]) * 1.0 / ((deltax)**2) + (u_old[i][j+1]-2*u_old[i][j]+u_old[i][j-1]) * 1.0 / ((deltay)**2)) for j in range(1,n)] for i in range(1,ms)])
    DIFV[1:ms, 1:n] = np.array([[nu * ((v_old[i+1][j]-2*v_old[i][j]+v_old[i-1][j]) * 1.0 / ((deltax)**2) + (v_old[i][j+1]-2*v_old[i][j]+v_old[i][j-1]) * 1.0 / ((deltay)**2)) for j in range(1,n)] for i in range(1,ms)])
    #u_old,v_old仮値設定②ナビエストークス方程式を解く
    u_old[1:ms-1, 1:n] = np.array([[u_old[i][j] + deltaT * (-(1.0/rho)*(p_old[i+1][j]-p_old[i][j])/deltax - CNVU[i][j] + DIFU[i][j]) for j in range(1, n)] for i in range(1, ms-1)])
    v_old[1:ms, 1:n-1] = np.array([[v_old[i][j] + deltaT * (-(1.0/rho)*(p_old[i][j+1]-p_old[i][j])/deltay - CNVV[i][j] + DIFV[i][j]) for j in range(1, n-1)] for i in range(1, ms)])
    #圧力補正ループ
    m1 = 1
    Dmax = M1 + 1
    while Dmax > M1 and np.min(deltap) > 1.0e-26 or np.max(deltap) == 0.0:
        #print "配列DIV内の最大値DmaxがMより大きい場合ループに入る"
        print str(value[1])
        print "t = " + str(t)
        print "m1 = " + str(m1)
        deltap[1:ms, 1:n] = np.array([[- rho * 1.0 / (2*deltaT) * (deltax*deltay) / (deltax**2+deltay**2) * (deltay * (u_old[i][j]-u_old[i-1][j]) + deltax * (v_old[i][j]-v_old[i][j-1])) for j in range(1, n)] for i in range(1, ms)])
        p_old[1:ms, 1:n] = np.array([[p_old[i][j] + 1.0 * omega * deltap[i][j] for j in range(1,n)] for i in range(1,ms)])
        u_old[1:ms-1, 1:n] = np.array([[u_old[i][j] + 1.0 * omega * (1.0/rho) * (deltaT/deltax) * deltap[i][j] for j in range(1,n)] for i in range(1,ms-1)])
        v_old[1:ms, 1:n-1] = np.array([[v_old[i][j] + 1.0 * omega * (1.0/rho) * (deltaT/deltay) * deltap[i][j] for j in range(1,n-1)] for i in range(1,ms)])
        u_old[1:ms-1, 1:n] = np.array([[u_old[i][j] - 1.0 * omega * (1.0/rho) * (deltaT/deltax) * deltap[i+1][j] for j in range(1,n)] for i in range(1,ms-1)])
        v_old[1:ms, 1:n-1] = np.array([[v_old[i][j] - 1.0 * omega * (1.0/rho) * (deltaT/deltay) * deltap[i][j+1] for j in range(1,n-1)] for i in range(1,ms)])
        boundary_condition()
        #連続の式の収束条件
        DIV[1:ms,1:n] = np.array([[abs((u_old[i][j] - u_old[i-1][j])*1.0/deltax + (v_old[i][j] - v_old[i][j-1])*1.0/deltay) for j in range(1,n)] for i in range(1,ms)])
        #---通常モード---
        Dmax = np.max(DIV)
        print "Dmax = " + str(Dmax)
        print "deltap = " + str(np.max(deltap))
        process_time = time.time() - start_time
        print "process_time = " + str(process_time)
        if m1 % 10000 == 0:
            csvout()
            graph()
        m1 += 1
        Dmax = 0#強制的ループ終了用

    #csvファイルで出力
    if t == deltaT or int(t/deltaT) % 25 == 0:
        csvout()
        graph()
    if t == deltaT or int(t/deltaT) % 50 == 0:
        if m1 % 10000 == 0:
            LineMessage()
            LineFigure("velocity(t=" + str(t) +",m1="+str(m1)+ ").png")
            LineFigure("pressure(t=" + str(t) +",m1="+str(m1)+ ").png")
        else:
            LineMessage()
            LineFigure("velocity(t=" + str(t) + ").png")
            LineFigure("pressure(t=" + str(t) + ").png")

    #時間を進める
    t = t + deltaT
print "Calculation ends"
LineMessage()
LineFigure("velocity(t=" + str(t) + ").png")
Linefigure("pressure(t=" + str(t) + ").png")
