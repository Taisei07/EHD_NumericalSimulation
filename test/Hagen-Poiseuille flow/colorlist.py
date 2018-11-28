# -*- coding: utf-8 -*-
    i = 0
    j = 0
    while 0 <= i <= n:
        while 0 <= j <= ms:
            velocity_out[i][j] = math.sqrt((u_out[i][j])**2+(v_out[i][j])**2)
            j += 1
        j = 0
        i += 1
    #色付け配列を設定
    i = 0
    j = 0
    while 0 <= i <= n:
        while 0 <= j <= ms:
            if velocity_out[i][j] <= (np.max(velocity_out)-np.min(velocity_out))/20*1:
                C_out[i][j] = 'b'
            elif (np.max(velocity_out)-np.min(velocity_out))/20*1 < velocity_out[i][j] <= (np.max(velocity_out)-np.min(velocity_out))/20*2:
                C_out[i][j] = 'navy'
            elif (np.max(velocity_out)-np.min(velocity_out))/20*2 < velocity_out[i][j] <= (np.max(velocity_out)-np.min(velocity_out))/20*3:
                C_out[i][j] = 'mediumblue'
            elif (np.max(velocity_out)-np.min(velocity_out))/20*3 < velocity_out[i][j] <= (np.max(velocity_out)-np.min(velocity_out))/20*4:
                C_out[i][j] = 'royalblue'
            elif (np.max(velocity_out)-np.min(velocity_out))/20*4 < velocity_out[i][j] <= (np.max(velocity_out)-np.min(velocity_out))/20*5:
                C_out[i][j] = 'steelblue'
            elif (np.max(velocity_out)-np.min(velocity_out))/20*5 < velocity_out[i][j] <= (np.max(velocity_out)-np.min(velocity_out))/20*6:
                C_out[i][j] = 'deepskyblue'
            elif (np.max(velocity_out)-np.min(velocity_out))/20*6 < velocity_out[i][j] <= (np.max(velocity_out)-np.min(velocity_out))/20*7:
                C_out[i][j] = 'darkturquoise'
            elif (np.max(velocity_out)-np.min(velocity_out))/20*7 < velocity_out[i][j] <= (np.max(velocity_out)-np.min(velocity_out))/20*8:
                C_out[i][j] = 'darkcyan'
            elif (np.max(velocity_out)-np.min(velocity_out))/20*8 < velocity_out[i][j] <= (np.max(velocity_out)-np.min(velocity_out))/20*9:
                C_out[i][j] = 'seagreen'
            elif (np.max(velocity_out)-np.min(velocity_out))/20*9 < velocity_out[i][j] <= (np.max(velocity_out)-np.min(velocity_out))/20*10:
                C_out[i][j] = 'limegreen'
            elif (np.max(velocity_out)-np.min(velocity_out))/20*10 < velocity_out[i][j] <= (np.max(velocity_out)-np.min(velocity_out))/20*11:
                C_out[i][j] = 'olivedrab'
            elif (np.max(velocity_out)-np.min(velocity_out))/20*11 < velocity_out[i][j] <= (np.max(velocity_out)-np.min(velocity_out))/20*12:
                C_out[i][j] = 'yellowgreen'
            elif (np.max(velocity_out)-np.min(velocity_out))/20*12 < velocity_out[i][j] <= (np.max(velocity_out)-np.min(velocity_out))/20*13:
                C_out[i][j] = 'yellow'
            elif (np.max(velocity_out)-np.min(velocity_out))/20*13 < velocity_out[i][j] <= (np.max(velocity_out)-np.min(velocity_out))/20*14:
                C_out[i][j] = 'khaki'
            elif (np.max(velocity_out)-np.min(velocity_out))/20*14 < velocity_out[i][j] <= (np.max(velocity_out)-np.min(velocity_out))/20*15:
                C_out[i][j] = 'gold'
            elif (np.max(velocity_out)-np.min(velocity_out))/20*15 < velocity_out[i][j] <= (np.max(velocity_out)-np.min(velocity_out))/20*16:
                C_out[i][j] = 'sandybrown'
            elif (np.max(velocity_out)-np.min(velocity_out))/20*16 < velocity_out[i][j] <= (np.max(velocity_out)-np.min(velocity_out))/20*17:
                C_out[i][j] = 'darkorange'
            elif (np.max(velocity_out)-np.min(velocity_out))/20*17 < velocity_out[i][j] <= (np.max(velocity_out)-np.min(velocity_out))/20*18:
                C_out[i][j] = 'tomato'
            elif (np.max(velocity_out)-np.min(velocity_out))/20*18 < velocity_out[i][j] <= (np.max(velocity_out)-np.min(velocity_out))/20*19:
                C_out[i][j] = 'orangered'
            elif (np.max(velocity_out)-np.min(velocity_out))/20*19 < velocity_out[i][j] <= (np.max(velocity_out)-np.min(velocity_out))/20*20:
                C_out[i][j] = 'red'
            j += 1
        j = 0
        i += 1
