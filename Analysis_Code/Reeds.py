import numpy as np

# 初始化变量
N = 14480  # 原子数量
Nc = 5  # number of correlation steps (a larger number gives a finer resolution)
Nf = Nc * 10  # The calculation numbers
fileName = "-OH20.txt"

# 打开轨迹文件
with open(fileName, "r") as f:
    # 读取轨迹文件
    p_all = np.zeros((Nf, N, 3))  # 存储所有原子的坐标
    for i in range(Nf):
        line = f.readline()
        if line.startswith("ITEM: ATOMS"):
            for k in range(N):
                line = f.readline().split()
                p_all[i, k] = [float(line[3]), float(line[4]), float(line[5])]

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

# 计算旋转半径均方
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

# 计算链端到端距离归一化均方
Re2x, Re2y, Re2z = np.zeros(Nf), np.zeros(Nf), np.zeros(Nf)
for i in range(Nf):
    for k in range(N - 1):
        dx, dy, dz = mic(p_all[i, k + 1, 0] - p_all[i, k, 0], p_all[i, k + 1, 1] - p_all[i, k, 1], p_all[i, k + 1, 2] - p_all[i, k, 2])
        Re2x[i] += dx ** 2
        Re2y[i] += dy ** 2
        Re2z[i] += dz ** 2

# 计算总的链端到端距离归一化均方
total_Reed2 = (Re2x + Re2y + Re2z) / (N - 1)

# 保存结果
np.savetxt("rg2.txt", np.array([Rg2x, Rg2y, Rg2z]).T, fmt="%.5f")
np.savetxt("re2.txt", np.array([Re2x, Re2y, Re2z]).T, fmt="%.5f")