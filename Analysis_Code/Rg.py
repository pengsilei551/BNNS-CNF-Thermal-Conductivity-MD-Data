import numpy as np
import matplotlib
matplotlib.use('TKAgg')
import matplotlib.pyplot as plt
import time

# 初始化变量
N = 16680  # number of atoms
Nc = 5  # number of correlation steps (a larger number gives a finer resolution)
dt = 0.0005  # time interval in units of ps (its inverse is roughly the maximum frequency attainable)
omega = np.arange(1, 380.5, 0.5)  # angular frequency in units of THz
nu = omega / 2 / np.pi  # omega = 2 * pi * nu, while nu is the frequency range
Nf = Nc * 10  # The calculation numbers
num_frame = Nf  # Frame numbers equals to calculation numbers
fileName = "-OH40.txt"  # Read the dump file

# 打开轨迹文件
with open(fileName, "r") as f:
    # 读取轨迹文件
    p_all = np.zeros((Nf, N, 3))  # 存储所有原子的坐标
    for i in range(Nf):
        line = f.readline()
        if line.startswith("ITEM: ATOMS"):
            for k in range(N):
                line = f.readline().split()
                if len(line) >= 5:  # 确保有足够的字段
                    # 坐标是第3、4、5个字段（从0开始计数）
                    p_all[i, k] = [float(line[2]), float(line[3]), float(line[4])]
                else:
                    print(f"Warning: Line {f.tell()} has insufficient data. Skipping this line.")
                    # 如果数据不足，用0填充
                    p_all[i, k] = [0.0, 0.0, 0.0]

# 计算质心
def mic(dx, dy, dz):
    Lx = Ly = Lz = 45  # 盒子尺寸
    dx = dx - Lx * np.round(dx / Lx)
    dy = dy - Ly * np.round(dy / Ly)
    dz = dz - Lz * np.round(dz / Lz)
    return dx, dy, dz

Rcmx, Rcmy, Rcmz = np.zeros(Nf), np.zeros(Nf), np.zeros(Nf)
for i in range(Nf):
    for k in range(N):
        Rcmx[i] += p_all[i, k, 0]
        Rcmy[i] += p_all[i, k, 1]
        Rcmz[i] += p_all[i, k, 2]
    Rcmx[i] /= N
    Rcmy[i] /= N
    Rcmz[i] /= N

# 计算均方回转半径（Rgs）
Rg2x, Rg2y, Rg2z = np.zeros(Nf), np.zeros(Nf), np.zeros(Nf)
for i in range(Nf):
    for k in range(N):
        dx, dy, dz = mic(p_all[i, k, 0] - Rcmx[i], p_all[i, k, 1] - Rcmy[i], p_all[i, k, 2] - Rcmz[i])
        Rg2x[i] += dx ** 2
        Rg2y[i] += dy ** 2
        Rg2z[i] += dz ** 2
    Rg2x[i] /= N
    Rg2y[i] /= N
    Rg2z[i] /= N

# 计算均方末端距离（Reeds）
Re2x, Re2y, Re2z = np.zeros(Nf), np.zeros(Nf), np.zeros(Nf)
for i in range(Nf):
    for k in range(N - 1):
        dx, dy, dz = mic(p_all[i, k + 1, 0] - p_all[i, k, 0], p_all[i, k + 1, 1] - p_all[i, k, 1], p_all[i, k + 1, 2] - p_all[i, k, 2])
        Re2x[i] += dx ** 2
        Re2y[i] += dy ** 2
        Re2z[i] += dz ** 2

# 保存结果
np.savetxt("rg2.txt", np.array([Rg2x, Rg2y, Rg2z]).T, fmt="%.5f")
np.savetxt("re2.txt", np.array([Re2x, Re2y, Re2z]).T, fmt="%.5f")

# 计算总的均方回转半径
total_Rg2 = (Rg2x + Rg2y + Rg2z) / N
average_Rg2 = total_Rg2.mean()

print("Average Rg2 over all frames:", average_Rg2)

# 计算总的均方末端距离
total_Reed2 = (Re2x + Re2y + Re2z) / (N - 1)

# 计算所有帧的平均均方末端距离
average_Reed2 = np.mean(total_Reed2)

print("Average Reed2 over all frames:", average_Reed2)