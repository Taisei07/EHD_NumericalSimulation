# -*- coding: utf-8 -*-

import math

import numpy as np
import itertools
from matplotlib import pyplot as plt
from matplotlib import animation


def update(frame, x, y):
    """グラフを更新するための関数"""
    # 現在のグラフを消去する
    plt.cla()
    # 折れ線グラフを再描画する
    plt.plot(x, y)
    plt.cla()


def main():
    # 描画領域
    fig = plt.figure(figsize=(10, 6))
    # 描画するデータ (最初は空っぽ)
    x = []
    y = []

    anime = animation.FuncAnimation(fig, update, fargs = (x, y), interval = 10, frames = itertools.count(0, 0.1))

    # グラフを表示する
    plt.show()


if __name__ == '__main__':
    main()
