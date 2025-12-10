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

# 计算平面法线
def calculate_normal_vector(v1, v2):
    return np.cross(v1, v2)

# 主函数
def main():
    dump_file = 'unshear33.8.txt'  # 请替换为您的dump文件路径
    atom_data = read_dump(dump_file)

    # 找到所有的B原子和N原子
    b_atoms = [atom for atom in atom_data if atom['type'] == 1]
    n_atoms = [atom for atom in atom_data if atom['type'] == 7]

    # 选择三个不共线的原子来定义平面
    if len(b_atoms) >= 3:
        # 选择前三个B原子来定义平面
        a1 = np.array([b_atoms[0]['x'], b_atoms[0]['y'], b_atoms[0]['z']])
        a2 = np.array([b_atoms[1]['x'], b_atoms[1]['y'], b_atoms[1]['z']])
        a3 = np.array([b_atoms[2]['x'], b_atoms[2]['y'], b_atoms[2]['z']])

        # 计算两个向量
        v1 = a2 - a1
        v2 = a3 - a1

        # 计算法线向量
        normal_vector = calculate_normal_vector(v1, v2)

        # 计算法线与X方向的夹角
        reference_vector = np.array([1, 0, 0])
        angle = calculate_angle(normal_vector, reference_vector)
        print(f"Orientation angle of h-BN sheet with X-axis: {angle} degrees")

if __name__ == "__main__":
    main()