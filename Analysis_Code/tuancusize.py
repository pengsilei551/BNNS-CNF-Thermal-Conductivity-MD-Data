def read_lammps_output(file_path):
    cluster_ids = {}
    with open(file_path, 'r') as file:
        for line in file:
            if 'ITEM: ATOMS' in line:
                break
        # 跳过头部信息，直到找到原子数据的开始
        for line in file:
            parts = line.split()
            if len(parts) < 4:
                continue  # 跳过不完整的行
            try:
                atom_id = int(parts[0])
                cluster_id = int(parts[1])
                cluster_ids[cluster_id] = cluster_ids.get(cluster_id, 0) + 1
            except ValueError:
                continue  # 如果转换失败，跳过这行
    return cluster_ids


def calculate_relative_size(cluster_ids):
    if not cluster_ids:
        return 0
    total_atoms = sum(cluster_ids.values())
    max_cluster_size = max(cluster_ids.values()) if cluster_ids else 0
    relative_size = max_cluster_size / total_atoms if total_atoms > 0 else 0
    return relative_size


def main():
    file_path = 'shear'  # 替换为你的文件路径
    cluster_ids = read_lammps_output(file_path)
    relative_size = calculate_relative_size(cluster_ids)

    print(f"Total number of clusters: {len(cluster_ids)}")
    for cluster_id, size in cluster_ids.items():
        print(f"Cluster {cluster_id}: {size} atoms")
    print(f"Relative size of the largest cluster: {relative_size:.4f}")


if __name__ == "__main__":
    main()