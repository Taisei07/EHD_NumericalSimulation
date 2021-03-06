# -*- coding: utf-8 -*-

import numpy as np
from numba import jit
import math
import matplotlib.pyplot as plt
import matplotlib.patches as patches
import json
import requests
import os
import csv
import sys
value = sys.argv
os.chdir('/media/pascal/HD-GDU3/Tajima_backup/EHD/result')
os.mkdir(str(value[1]))
import time
start_time = time.time()
def time_count():
    process_time = time.time() - start_time
    print "process_time = " + str(process_time)
import requests

#slack通知
def slack_mention():
    process_time = time.time() - start_time
    requests.post('https://hooks.slack.com/services/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', data = json.dumps({
        'text': str(value[1]) + '\n' + 'process_time : ' + str(process_time) + '\n' + 't : ' + str(t) , # 投稿するテキスト
        'username': u'ghost', # 投稿のユーザー名
        'icon_emoji': u':ghost:', # 投稿のプロフィール画像に入れる絵文字
        'link_names': 1, # メンションを有効にする
    }))

def figure_upload(A):
    os.chdir(str(value[1]))
    TOKEN = 'xoxb-335429545841-495712943572-ogOBXk5mgALgEKTaDkcick7e'
    CHANNEL = 'CELQHE10X'
    files = {'file': open(A, 'rb')}
    param = {
    'token':TOKEN,
    'channels':CHANNEL,
    'filename':A,
    'initial_comment': str(value[1]) + '_' + A,
    'title': A
    }
    requests.post(url="https://slack.com/api/files.upload",params=param, files=files)
    os.chdir('../')

#物性値
print "物性値入力"
rho = 950#input("rho(流体密度)[kg/m^3] = ")
nu = 0.000096#input("nu（動粘性係数）[m^2/s] = ")
epsilon = 0.000000000099#input("epsilon(流体の誘電率)[F/m] = ")
K = 0.000000001#input("K(イオン移動度)[m^2/Vs] = ")
Di = 0#input("Di(イオン拡散係数) = ")
sigma = 0.000000001#input("sigma(導電率)[A/Vm] = ")

#定数
print "定数入力"
H = 0.002#input("H(流れ場y方向長さ)[m] = ")
L = 0.006#input("L(流れ場x方向長さ)[m] = ")
T = input("T(移流時間)[s] = ")
deltaT = 0.00001#input("deltaT(時間刻み)[s] = ")
deltax = input("deltax(x方向要素間距離)[m] = ")
deltay = input("deltay(y方向要素間距離)[m] = ")
omega = 0.5#input("omega(緩和係数) = ")
M1 = 0.000001#input("M1(連続の式収束条件) = ")
M2 = 0.00000001#input("M2(電位phiの収束条件) = ")
electrode_number = input("electrode_number(電極の数,2or4or8)[個] = ")
E_x = 0.0019#input("E_x(electrode上の点Eのx座標)[m] = ")
F_x = 0.0029#input("F_x(electrode上の点Cのx座標)[m] = ")
G_x = 0.0031#input("G_x(electrode上の点Dのx座標)[m] = ")
H_x = 0.0041#input("H_x(electrode上の点Eのx座標)[m] = ")
if electrode_number == 4 or electrode_number == 8:
    electrode_pattern = input("electrode_pattern(電極配置のパターン,line or topandbottom) = ")
    if electrode_pattern == "line":
        L = 0.0084
        I_x = 0.0043#input("I_x(electrode上の点Iのx座標)[m] = ")
        J_x = 0.0053#input("J_x(electrode上の点Jのx座標)[m] = ")
        K_x = 0.0055#input("K_x(electrode上の点Kのx座標)[m] = ")
        L_x = 0.0065#input("L_x(electrode上の点Lのx座標)[m] = ")
phi_electrodeEF = 0#input("phi_electrodeEF(電極EFの電位)[V] = ")
phi_electrodeGH = input("phi_electrodeGH(電極GHの電位)[V] = ")
if electrode_number == 4 or electrode_number == 8:
    phi_electrodeIJ = 0#input("phi_electrodeIJ(電極IJの電位)[V] = ")
    phi_electrodeKL = input("phi_electrodeKL(電極KLの電位)[V] = ")


#物性値、定数の出力
constant_list = [\
["密度rho[kg/m^3]",rho],\
["動粘性係数nu[m^2/s]",nu],\
["流体の誘電率epsilon[F/m]",epsilon],\
["イオン移動度K[m^2/Vs]",K],\
["イオン拡散係数Di",Di],\
["導電率sigma[A/Vm]",sigma],\
["流れ場y方向長さH[m]",H],\
["流れ場x方向長さL[m]",L],\
["移流時間T[s]",T],\
["時間刻みdeltaT[s]",deltaT],\
["x方向要素間距離deltax[m]",deltax],\
["y方向要素間距離deltay[m]",deltay],\
["緩和係数omega",omega],\
["連続の式収束条件M1",M1],\
["電位phiの収束条件M2",M2],\
["電極EFの電位phi_electrodeEF[V]",phi_electrodeEF],\
["電極EHの電位phi_electrodeGH[V]",phi_electrodeGH]]
with open(os.path.join(str(value[1]),"constant.csv"), 'w') as file:
    writer = csv.writer(file, lineterminator = '\n')
    writer.writerows(constant_list)

#流れ場の座標(ms,n)
ms = int(L / deltax + 1)
n = int(H / deltay + 1)
print "ms :" + str(ms)
print "n :" + str(n)

#格子定義点(oldはtimestep=n, newはtimestep=n+1)　＊座標(ms,n)に1を足している
u_old = np.array([[0.0] * (n+1) for i in range(ms+1)])#ms:x方向,n:y方向
v_old = np.array([[0.0] * (n+1) for i in range(ms+1)])
p_old = np.array([[0.0] * (n+1) for i in range(ms+1)])
phi = np.array([[0.0] * (n+1) for i in range(ms+1)])
phi_old = np.array([[0.0] * (n+1) for i in range(ms+1)])
Ex = np.array([[0.0] * (n+1) for i in range(ms+1)])
Ey = np.array([[0.0] * (n+1) for i in range(ms+1)])
E = np.array([[0.0] * (n+1) for i in range(ms+1)])
q = np.array([[0.0] * (n+1) for i in range(ms+1)])
Fx = np.array([[0.0] * (n+1) for i in range(ms+1)])
Fy = np.array([[0.0] * (n+1) for i in range(ms+1)])
F = np.array([[0.0] * (n+1) for i in range(ms+1)])

#対流項CNVと粘性項DIF配列設定
CNVU = np.array([[0.0] * (n+1) for i in range(ms+1)])
CNVV = np.array([[0.0] * (n+1) for i in range(ms+1)])
DIFU = np.array([[0.0] * (n+1) for i in range(ms+1)])
DIFV = np.array([[0.0] * (n+1) for i in range(ms+1)])
#イオン移動度の項IMOB、イオン運動量の項IMOM、イオン拡散の項IDIF、導電率の項ICODの配列設定
IMOB = np.array([[0.0] * (n+1) for i in range(ms+1)])
IMOM = np.array([[0.0] * (n+1) for i in range(ms+1)])
IDIF = np.array([[0.0] * (n+1) for i in range(ms+1)])
ICOD = np.array([[0.0] * (n+1) for i in range(ms+1)])
#連続の式DIV・電位DPHの配列設定
DIV = np.array([[0.0] * (n+1) for i in range(ms+1)])

#境界条件の設定
def boundary_condition():
    #BoundaryAD
    for j in range(n):
        u_old[0][j] = u_old[1][j]
        v_old[0][j] = v_old[1][j]
        p_old[0][j] = p_old[1][j]
        j += 1
    #WallAB
    for i in range(ms):
        u_old[i][0] = -u_old[i][1]
        v_old[i][0] = 0.0
        p_old[i][0] = p_old[i][1]
        i += 1
    #WallCD
    for i in range(ms):
        u_old[i][n] = -u_old[i][n-1]
        v_old[i][n-1] = 0.0
        p_old[i][n] = p_old[i][n-1]
        i += 1
    #BoundaryBC
    for j in range(n):
        u_old[ms-1][j] = u_old[ms-2][j]
        v_old[ms][j] = v_old[ms-1][j]
        p_old[ms][j] = p_old[ms-1][j]
        j += 1

def boundary_condition_phi():
    #BoundaryAD
    for j in range(n):
        phi[0][j] = phi[1][j]
        j += 1
    #WallAB
    for i in range(ms):
        phi[i][0] = phi[i][1]
        i += 1
    #WallCD
    for i in range(ms):
        phi[i][n] = phi[i][n-1]
        i += 1
    #BoundaryBC
    for j in range(n):
        phi[ms][j] = phi[ms-1][j]
        j += 1
    #ElectrodeEF
    for i in range(int(E_x / deltax + 1), int(F_x / deltax + 1)+1):
        phi[i][0] = phi_electrodeEF
        i += 1
    #ElectrodeGH
    for i in range(int(G_x / deltax + 1), int(H_x / deltax + 1)+1):
        phi[i][0] = phi_electrodeGH
        i += 1

    if electrode_number == 4 or electrode_number == 8:
        if electrode_pattern == "line":
            #ElectrodeIJ
            for i in range(int(I_x / deltax + 1), int(J_x / deltax + 1)+1):
                phi[i][0] = phi_electrodeIJ
                i += 1
            #ElectrodeKL
            for i in range(int(K_x / deltax + 1), int(L_x / deltax + 1)+1):
                phi[i][0] = phi_electrodeKL
                i += 1
        elif electrode_pattern == "topandbottom":
            #ElectrodeIJ
            for i in range(int(E_x / deltax + 1), int(F_x / deltax + 1)+1):
                phi[i][n] = phi_electrodeIJ
                i += 1
            #ElectrodeKL
            for i in range(int(G_x / deltax + 1), int(H_x / deltax + 1)+1):
                phi[i][n] = phi_electrodeKL
                i += 1

def boundary_condition_electrode(A,B):
    for i in range(int(A / deltax + 1), int(B / deltax + 1)+1):
        q[i][1] = - epsilon * (phi[i][2]-phi[i][1]) / (deltay**2)
        q[i][0] = 0
        i += 1

def boundary_condition_q():
    #BoundaryAD
    for j in range(n):
        q[0][j] = q[0][j]
        j += 1
    #WallAB
    for i in range(ms):
        q[i][0] = q[i][1]
        i += 1
    #WallCD
    for i in range(ms):
        q[i][n] = q[i][n-1]
        i += 1
    #BoundaryBC
    for j in range(n):
        q[ms][j] = q[ms-1][j]
        j += 1
    #ElectrodeEF
    boundary_condition_electrode(E_x,F_x)
    #ElectrodeGH
    boundary_condition_electrode(G_x,H_x)

    if electrode_number == 4 or electrode_number == 8:
        if electrode_pattern == "line":
            #ElectrodeIJ
            boundary_condition_electrode(I_x,J_x)
            #ElectrodeKL
            boundary_condition_electrode(K_x,L_x)
        elif electrode_pattern == "topandbottom":
            #ElectrodeIJ
            for i in range(int(E_x / deltax + 1), int(F_x / deltax + 1)+1):
                q[i][n-1] = - epsilon * (phi[i][n-2]-phi[i][n-1]) / (deltay**2)
                q[i][n] = 0
                i += 1
            #ElectrodeKL
            for i in range(int(G_x / deltax + 1), int(H_x / deltax + 1)+1):
                q[i][n-1] = - epsilon * (phi[i][n-2]-phi[i][n-1]) / (deltay**2)
                q[i][n] = 0
                i += 1

#初期におけるDIV
@jit
def DIV_calculation():
    for i in range(1, ms):
        for j in range(1, n):
            DIV[i][j] = abs((u_old[i][j] - u_old[i-1][j])*1.0/deltax + (v_old[i][j] - v_old[i][j-1])*1.0/deltay)
            j += 1
        j = 1
        i += 1

#過去の電位phiを更新する
@jit
def phi_old_def():
    for i in range(ms+1):
        for j in range(n+1):
            phi_old[i][j] = phi[i][j]
            j += 1
        j = 0
        i += 1

#電場強度Eの計算
@jit
def E_calculation():
    for i in range(ms):
        for j in range(n):
            Ex[i][j] = -(phi[i+1][j]-phi[i][j])/deltax
            Ey[i][j] = -(phi[i][j+1]-phi[i][j])/deltay
            E[i][j] = np.sqrt(Ex[i][j]**2+Ey[i][j]**2)
            j += 1
        j = 0
        i += 1

#電位phiの計算
@jit
def phi_calculation():
    DPmax = M2 + 1
    m1 = 1
    while DPmax > M2:
        print str(value[1])
        print "t = " + str(t)
        print "m1 = " + str(m1)
        phi_old_def()
        for j in xrange(1, n):
            for i in xrange(1, ms):
                phi[i][j] = 1.0 * (deltax * deltay)**2 / (2 * ((deltax**2)+(deltay**2))) * (1.0 * q[i][j] / epsilon + (phi[i+1][j]+phi[i-1][j])/(deltax**2) + (phi[i][j+1]+phi[i][j-1])/(deltay**2))
                i += 1
            i = 1
            j += 1
        boundary_condition_phi()
        DPmax = np.max(phi-phi_old)
        print "DPmax = " + str(DPmax)
        m1 += 1
    E_calculation()

#電気的な力Fの計算
@jit
def F_calculation():
    for i in range(ms+1):
        for j in range(n+1):
            F[i][j] = np.sqrt(Fx[i][j]**2+Fy[i][j]**2)
            j += 1
        j = 1
        i += 1

#電荷保存則の計算
@jit
def q_calculation():
    for j in range(1, n):
        for i in range(1, ms):
            IMOB[i][j] = 1.0 * q[i][j] * ((Ex[i][j]-Ex[i-1][j])/deltax + (Ey[i][j]-Ey[i][j-1])/deltay) + 1.0 * ((Ex[i+1][j]+Ex[i][j])/2*(q[i+1][j]-q[i][j])/2/deltax + (Ey[i][j+1]+Ey[i][j])/2*(q[i][j+1]-q[i][j])/2/deltay)
            IMOM[i][j] = 1.0 * q[i][j] * ((u_old[i][j] - u_old[i-1][j]) / deltax + 1.0 * (v_old[i][j] - v_old[i][j-1]) / deltay) + 1.0 * (u_old[i][j] + u_old[i-1][j]) / 2 * (q[i+1][j] - q[i-1][j]) / (2*deltax) + 1.0 * (v_old[i][j] + v_old[i][j-1]) / 2 * (q[i][j+1]-q[i][j-1]) / (2*deltay)
            IDIF[i][j] = 1.0 * (q[i+1][j] - 2 * q[i][j] + q[i-1][j]) / (deltax**2) + 1.0 * (q[i][j+1] - 2 * q[i][j] + q[i][j-1]) / (deltay**2)
            ICOD[i][j] = 1.0 * ((Ex[i][j]-Ex[i-1][j])/deltax + (Ey[i][j]-Ey[i][j-1])/deltay)
            i += 1
        i = 1
        j += 1
    for j in range(1, n):
        for i in range(1, ms):
            q[i][j] = q[i][j] - deltaT * (K * IMOB[i][j] + IMOM[i][j] - Di * IDIF[i][j] + sigma * ICOD[i][j])
            i += 1
        i = 1
        j += 1

#電気的な力FX,Fyの配列設定
@jit
def FxFy_calculation():
    for i in range(1, ms):
        for j in range(1, n):
            Fx[i][j] = - 1.0 * (q[i+1][j]+q[i][j]) / 2.0 * Ex[i][j]
            Fy[i][j] = - 1.0 * (q[i][j+1]+q[i][j]) / 2.0 * Ey[i][j]
            j += 1
        j = 1
        i += 1

#u_old,v_old仮値設定①粘性項・対流項配列の設定
@jit
def temporary_velocity():
    for i in range(1, ms):
        for j in range(1, n):
            CNVU[i][j] = (u_old[i][j] * (u_old[i+1][j]-u_old[i-1][j]) * 1.0 /(2*deltax) - abs(u_old[i][j]) * (u_old[i+1][j]-2*u_old[i][j]+u_old[i-1][j]) * 1.0 / (2*deltax)) - (v_old[i][j]*(u_old[i][j+1] - u_old[i][j-1]) * 1.0 / (2 * deltay) - abs(v_old[i][j])*(u_old[i][j+1]-2*u_old[i][j]+u_old[i][j-1]) * 1.0 / (2 * deltay))
            CNVV[i][j] = (u_old[i][j] * (v_old[i+1][j]-v_old[i-1][j]) * 1.0 /(2*deltax) - abs(u_old[i][j]) * (v_old[i+1][j]-2*v_old[i][j]+v_old[i-1][j]) * 1.0 / (2*deltax)) - (v_old[i][j]*(v_old[i][j+1] - v_old[i][j-1]) * 1.0 / (2 * deltay) - abs(v_old[i][j])*(v_old[i][j+1]-2*v_old[i][j]+v_old[i][j-1]) * 1.0 / (2 * deltay))
            DIFU[i][j] = nu * ((u_old[i+1][j]-2*u_old[i][j]+u_old[i-1][j]) * 1.0 / ((deltax)**2) + (u_old[i][j+1]-2*u_old[i][j]+u_old[i][j-1]) * 1.0 / ((deltay)**2))
            DIFV[i][j] = nu * ((v_old[i+1][j]-2*v_old[i][j]+v_old[i-1][j]) * 1.0 / ((deltax)**2) + (v_old[i][j+1]-2*v_old[i][j]+v_old[i][j-1]) * 1.0 / ((deltay)**2))
            j += 1
        j = 1
        i += 1

#u_old,v_old仮値設定②ナビエストークス方程式を解く
@jit
def navie_calculation():
    for i in range(1, ms-1):
        for j in range(1, n):
            u_old[i][j] = u_old[i][j] + deltaT * (-(1.0/rho)*(p_old[i+1][j]-p_old[i][j])/deltax - CNVU[i][j] + DIFU[i][j] + 1.0 / rho * Fx[i][j])
            j += 1
        j = 1
        i += 1
    for i in range(1, ms):
        for j in range(1, n-1):
            v_old[i][j] = v_old[i][j] + deltaT * (-(1.0/rho)*(p_old[i][j+1]-p_old[i][j])/deltay - CNVV[i][j] + DIFV[i][j] + 1.0 / rho * Fy[i][j])
            j += 1
        j = 1
        i += 1

#速度補正ループ
@jit
def velocity_correct():
    for i in range(1, ms):
        for j in range(1, n):
            deltap = - rho * 1.0 / (2*deltaT) * (deltax*deltay) / (deltax**2+deltay**2) * (deltay * (u_old[i][j]-u_old[i-1][j]) + deltax * (v_old[i][j]-v_old[i][j-1]))
            p_old[i][j] = p_old[i][j] + 1.0 * omega * deltap
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

#csvファイルで出力
def csvout01():
    phi_out = phi.transpose()
    q_out = q.transpose()
    E_out = E.transpose()
    #電気的な力Fの計算
    F_calculation()
    F_out = F.transpose()
    import csv
    with open(os.path.join(str(value[1]),"phi_(t="+str(t)+")"+".csv"), 'w') as file:
        writer = csv.writer(file, lineterminator = '\n')
        writer.writerows(phi_out)
    with open(os.path.join(str(value[1]),"q_(t="+str(t)+")"+".csv"), 'w') as file:
        writer = csv.writer(file, lineterminator = '\n')
        writer.writerows(q_out)
    with open(os.path.join(str(value[1]),"F(t="+str(t)+")"+".csv"), 'w') as file:
        writer = csv.writer(file, lineterminator = '\n')
        writer.writerows(F_out)
    with open(os.path.join(str(value[1]),"E(t="+str(t)+")"+".csv"), 'w') as file:
        writer = csv.writer(file, lineterminator = '\n')
        writer.writerows(E_out)

def csvout02():
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


#グラフを作成して保存する
#分布図用座標
X1 = np.array([[0.0] * (n+1) for i in range(ms+1)])
for j in range(n+1):
    for i in range(ms+1):
        X1[i][j] = deltax * (i-0.5)
        i += 1
    i = 0
    j += 1
Y1 = np.array([[0.0] * (n+1) for i in range(ms+1)])
for j in range(n+1):
    for i in range(ms+1):
        Y1[i][j] = deltay * (j-0.5)
        i += 1
    i = 0
    j += 1
X1_out = X1.transpose()
Y1_out = Y1.transpose()
#ベクトル用座標
X2 = np.array([[0.0] * (n+1) for i in range(ms+1)])
for j in range(n+1):
    for i in range(ms+1):
        X2[i][j] = deltax * i
        i += 1
    i = 0
    j += 1
Y2 = np.array([[0.0] * (n+1) for i in range(ms+1)])
for j in range(n+1):
    for i in range(ms+1):
        Y2[i][j] = deltay * j
        i += 1
    i = 0
    j += 1
X2_out = X2.transpose()
Y2_out = Y2.transpose()

def fig_electrode():
    ax = plt.axes()
    e1 = patches.Rectangle(xy=(E_x, -0.00005), width=F_x-E_x, height=0.00005, fc='y')
    ax.add_patch(e1)
    e2 = patches.Rectangle(xy=(G_x, -0.00005), width=H_x-G_x, height=0.00005, fc='y')
    ax.add_patch(e2)
    r = patches.Rectangle(xy=(0, 0), width=L, height=H, ec='#000000', fill=False)
    ax.add_patch(r)
    if electrode_number == 4 or electrode_number == 8:
        if electrode_pattern == "line":
            e3 = patches.Rectangle(xy=(I_x, -0.00005), width=J_x-I_x, height=0.00005, fc='y')
            ax.add_patch(e3)
            e4 = patches.Rectangle(xy=(K_x, -0.00005), width=L_x-K_x, height=0.00005, fc='y')
            ax.add_patch(e4)
        elif electrode_pattern == "topandbottom":
            e3 = patches.Rectangle(xy=(E_x, H), width=F_x-E_x, height=0.00005, fc='y')
            ax.add_patch(e3)
            e4 = patches.Rectangle(xy=(G_x, H), width=H_x-G_x, height=0.00005, fc='y')
            ax.add_patch(e4)

def graph_enlarge():
    if electrode_number == 2:
        plt.xlim(-1.0*L/20 + E_x, H_x + 1.0*L/20)
        plt.ylim(-1.0*H/20, 1.0*H/4)
    elif electrode_number == 4:
        if electrode_pattern == 'line':
            plt.xlim(-1.0*L/20 + E_x, L_x + 1.0*L/20)
            plt.ylim(-1.0*H/20, 1.0*H/4)
        elif electrode_pattern == 'topandbottom':
            plt.xlim(-1.0*L/20 + E_x, H_x + 1.0*L/20)
            plt.ylim(-1.0*H/20, 11.0*H/10)
def graph01():
    #配列変換・設定
    phi_out = phi.transpose()
    q_out = q.transpose()
    E_out = E.transpose()
    Ex_out = Ex.transpose()
    Ey_out = Ey.transpose()
    Fx_out = Fx.transpose()
    Fy_out = Fy.transpose()
    F_calculation()
    F_out = F.transpose()
    #ディレクトリ移動
    os.chdir(str(value[1]))
    #電荷密度分布作成
    plt.pcolor(X1_out, Y1_out, q_out)
    plt.colorbar()
    plt.axis('equal')
    plt.title('electricalcharge_distribution(t='+str(t)+')')
    plt.xlabel('x')
    plt.ylabel('y')
    plt.xlim(-1.0*L/10, 11.0*L/10)
    plt.ylim(-1.0*H/10, 11.0*H/10)
    plt.grid()
    fig_electrode()
    plt.savefig("electricalcharge(t=" + str(t) + ").png", dpi=600)
    plt.cla()
    plt.clf()
    plt.close()
    #電荷密度分布(拡大図)
    plt.pcolor(X1_out, Y1_out, q_out)
    plt.colorbar()
    plt.axis('equal')
    plt.title('electricalcharge_distribution(t='+str(t)+')')
    plt.xlabel('x')
    plt.ylabel('y')
    graph_enlarge()
    plt.grid()
    fig_electrode()
    plt.savefig("electricalcharge(enlarge, t=" + str(t) + ").png", dpi=600)
    plt.cla()
    plt.clf()
    plt.close()
    #電位分布作成
    plt.pcolor(X1_out, Y1_out, phi_out)
    plt.colorbar()
    plt.axis('equal')
    plt.title('electricalpotential_distribution(t='+str(t)+')')
    plt.xlabel('x')
    plt.ylabel('y')
    plt.xlim(-1.0*L/10, 11.0*L/10)
    plt.ylim(-1.0*H/10, 11.0*H/10)
    plt.grid()
    fig_electrode()
    plt.savefig("electricalpotential(t=" + str(t) + ").png", dpi=600)
    plt.cla()
    plt.clf()
    plt.close()
    #電場強度ベクトル作成
    plt.pcolor(X1_out, Y1_out, E_out)
    plt.colorbar()
    plt.quiver(X2_out, Y2_out, Ex_out, Ey_out, angles='xy', scale_units='xy', scale=np.max(E_out)*(1.0/(L*0.03/2)), headwidth=5, headlength=8, headaxislength=4)
    plt.axis('equal')
    plt.title('electrofield_vector(t='+str(t)+')')
    plt.xlabel('x')
    plt.ylabel('y')
    plt.xlim(-1.0*L/10, 11.0*L/10)
    plt.ylim(-1.0*H/10, 11.0*H/10)
    plt.grid()
    plt.draw()
    fig_electrode()
    plt.savefig("electrofield(t=" + str(t) + ").png", dpi=600)
    plt.cla()
    plt.clf()
    plt.close()
    #電気的な力Fベクトル作成
    plt.pcolor(X1_out, Y1_out, F_out)
    plt.colorbar()
    plt.quiver(X2_out, Y2_out, Fx_out, Fy_out, angles='xy', scale_units='xy', scale=np.max(F_out)*(1.0/(L*0.03/2)), headwidth=5, headlength=8, headaxislength=4)
    plt.axis('equal')
    plt.title('F_vector(t='+str(t)+')')
    plt.xlabel('x')
    plt.ylabel('y')
    plt.xlim(-1.0*L/10, 11.0*L/10)
    plt.ylim(-1.0*H/10, 11.0*H/10)
    plt.grid()
    plt.draw()
    fig_electrode()
    plt.savefig("F(t=" + str(t) + ").png", dpi=600)
    plt.cla()
    plt.clf()
    plt.close()
    #電気的な力Fベクトル(拡大図)
    plt.pcolor(X1_out, Y1_out, F_out)
    plt.colorbar()
    plt.quiver(X2_out, Y2_out, Fx_out, Fy_out, angles='xy', scale_units='xy', scale=np.max(F_out)*(1.0/(L*0.03/2)), headwidth=5, headlength=8, headaxislength=4)
    plt.axis('equal')
    plt.title('F_vector(t='+str(t)+')')
    plt.xlabel('x')
    plt.ylabel('y')
    graph_enlarge()
    plt.grid()
    plt.draw()
    fig_electrode()
    plt.savefig("F(enlarge, t=" + str(t) + ").png", dpi=600)
    plt.cla()
    plt.clf()
    plt.close()
    os.chdir('../')

def graph02():
    #配列変換・設定
    u_out = u_old.transpose()
    v_out = v_old.transpose()
    p_out = p_old.transpose()
    velocity_out = np.sqrt(u_out**2+v_out**2)
    #ディレクトリ移動
    os.chdir(str(value[1]))
    #速度ベクトル作成
    plt.pcolor(X1_out, Y1_out, velocity_out)
    plt.colorbar()
    plt.quiver(X2_out, Y2_out, u_out, v_out, angles='xy', scale_units='xy', scale=np.max(velocity_out)*(1.0/(L*0.03/2)), headwidth=5, headlength=8, headaxislength=4)
    plt.axis('equal')
    plt.title('velocity_vector(t='+str(t)+')')
    plt.xlabel('x')
    plt.ylabel('y')
    plt.xlim(-1.0*L/10, 11.0*L/10)
    plt.ylim(-1.0*H/10, 11.0*H/10)
    plt.grid()
    plt.draw()
    fig_electrode()
    plt.savefig("velocity(t=" + str(t) + ").png", dpi=600)
    plt.cla()
    plt.clf()
    plt.close()
    #速度ベクトル(拡大図)
    plt.pcolor(X1_out, Y1_out, velocity_out)
    plt.colorbar()
    plt.quiver(X2_out, Y2_out, u_out, v_out, angles='xy', scale_units='xy', scale=np.max(velocity_out)*(1.0/(L*0.03/2)), headwidth=5, headlength=8, headaxislength=4)
    plt.axis('equal')
    plt.title('velocity_vector(t='+str(t)+')')
    plt.xlabel('x')
    plt.ylabel('y')
    graph_enlarge()
    plt.grid()
    plt.draw()
    fig_electrode()
    plt.savefig("velocity(enlarge, t=" + str(t) + ").png", dpi=600)
    plt.cla()
    plt.clf()
    plt.close()
    #圧力分布作成
    plt.pcolor(X1_out, Y1_out, p_out)
    plt.colorbar()
    plt.axis('equal')
    plt.title('pressure_distribution(t='+str(t)+')')
    plt.xlabel('x')
    plt.ylabel('y')
    plt.xlim(-1.0*L/10, 11.0*L/10)
    plt.ylim(-1.0*H/10, 11.0*H/10)
    plt.grid()
    fig_electrode()
    plt.savefig("pressure(t=" + str(t) + ").png", dpi=600)
    plt.cla()
    plt.clf()
    plt.close()
    os.chdir('../')

def reset_uv():
    u_old = np.array([[0.0] * (n+1) for i in range(ms+1)])
    v_old = np.array([[0.0] * (n+1) for i in range(ms+1)])


print "Calculation starts"
#初期状態の計算
t = 0
m1 = 1
m2 = 1
boundary_condition_phi()
phi_calculation()
DIV_calculation()
csvout01()
csvout02()
graph01()
graph02()
slack_mention()
figure_upload("electricalpotential(t=" + str(t) + ").png")
figure_upload("electrofield(t=" + str(t) + ").png")
t = deltaT
for i in xrange(int(T/deltaT)):
    print "t =" + str(t)
    q_calculation()
    boundary_condition_q()
    phi_calculation()
    FxFy_calculation()
    csvout01()
    graph01()
    temporary_velocity()
    navie_calculation()

    #圧力補正ループ
    m2 = 1
    Dmax = M1 + 1
    while Dmax > M1 :
        #print "配列DIV内の最大値DmaxがMより大きい場合ループに入る"
        print str(value[1])
        print "t = " + str(t)
        print "m2 = " + str(m2)
        velocity_correct()
        boundary_condition()
        DIV_calculation()
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
        time_count()
        m2 += 1
        #Dmax = 0#強制的ループ終了用
    #csvファイルで出力
    if t == deltaT or int(t/deltaT) % 5 == 0:
        csvout02()
        graph02()
    if t == deltaT or int(t/deltaT) % 5 == 0:
        slack_mention()
        figure_upload("electricalcharge(t=" + str(t) + ").png")
        figure_upload("electricalcharge(enlarge, t=" + str(t) + ").png")
        figure_upload("electricalpotential(t=" + str(t) + ").png")
        figure_upload("electrofield(t=" + str(t) + ").png")
        figure_upload("F(t=" + str(t) + ").png")
        figure_upload("F(enlarge, t=" + str(t) + ").png")
        figure_upload("velocity(t=" + str(t) + ").png")
        figure_upload("velocity(enlarge, t=" + str(t) + ").png")
        figure_upload("pressure(t=" + str(t) + ").png")

    #時間を進める
    t = t + deltaT
print "Calculation ends"
slack_mention()
figure_upload("electricalcharge(t=" + str(t) + ").png")
figure_upload("electricalcharge(enlarge, t=" + str(t) + ").png")
figure_upload("electricalpotential(t=" + str(t) + ").png")
figure_upload("electrofield(t=" + str(t) + ").png")
figure_upload("F(t=" + str(t) + ").png")
figure_upload("F(enlarge, t=" + str(t) + ").png")
figure_upload("velocity(t=" + str(t) + ").png")
figure_upload("velocity(enlarge, t=" + str(t) + ").png")
figure_upload("pressure(t=" + str(t) + ").png")
