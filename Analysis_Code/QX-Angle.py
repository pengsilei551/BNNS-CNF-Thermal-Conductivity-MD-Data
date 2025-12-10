import numpy as np


# 读取LAMMPS dump文件
def read_dump(file_name):
    atom_data = []
    with open(file_name, 'r') as f:
        lines = f.readlines()

    # 跳过非ATOMS步骤
    in_atoms = False
    for line in lines:
        if line.startswith('ITEM: ATOMS'):
            in_atoms = True
            continue
        if in_atoms and line.strip() and not line.startswith('ITEM:'):
            parts = line.split()
            if len(parts) >= 5:
                id, type, x, y, z = int(parts[0]), int(parts[1]), float(parts[2]), float(parts[3]), float(parts[4])
                atom_data.append({'id': id, 'type': type, 'x': x, 'y': y, 'z': z})
        if line.startswith('ITEM: BOX BOUNDS'):
            in_atoms = False
    return atom_data


# 计算两个向量之间的夹角
def calculate_angle(v1, v2):
    dot_product = np.dot(v1, v2)
    magnitude_v1 = np.linalg.norm(v1)
    magnitude_v2 = np.linalg.norm(v2)
    angle = np.arccos(dot_product / (magnitude_v1 * magnitude_v2))
    return np.degrees(angle)


# 主函数
def main():
    dump_file = 'shear33.8.txt'  # 请替换为您的dump文件路径
    atom_data = read_dump(dump_file)

    # 找到所有的B原子和N原子
    b_atoms = [atom for atom in atom_data if atom['type'] == 1]  # 修正此处
    n_atoms = [atom for atom in atom_data if atom['type'] == 7]  # 修正此处

    # 计算B原子与N原子之间的取向角
    for b_atom in b_atoms:
        for n_atom in n_atoms:
            if b_atom['id'] != n_atom['id']:  # 确保不是同一个原子
                # 计算B和N原子之间的向量
                vector_bn = np.array([b_atom['x'] - n_atom['x'], b_atom['y'] - n_atom['y'], b_atom['z'] - n_atom['z']])
                # 计算取向角，这里只考虑X方向，所以向量是(1, 0, 0)
                angle = calculate_angle(vector_bn, (1, 0, 0))
                print(f"B atom {b_atom['id']} and N atom {n_atom['id']} angle with X-axis: {angle} degrees")


if __name__ == "__main__":
    main()