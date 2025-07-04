import numpy as np
import os
import ase.io # 确保 ase 库已安装并导入
import glob
import re

# 定义一个小的数值稳定性容差，用于接近零的体积
VOLUME_TOLERANCE = 1e-9

def calculate_cell_volume(lattice_str):
    """
    根据晶胞向量字符串计算晶胞体积。
    字符串格式通常是 "Ax Ay Az Bx By Bz Cx Cy Cz"。
    体积通过晶胞矩阵的行列式绝对值计算。
    """
    try:
        coords = np.array(list(map(float, lattice_str.split())))
        if len(coords) != 9:
            print(f"DEBUG(Volume): 晶胞字符串长度不正确 ({len(coords)}): '{lattice_str}'")
            return 0.0 # 视为无效晶胞，体积为零
        
        cell_matrix = coords.reshape(3, 3)
        volume = np.abs(np.linalg.det(cell_matrix))
        
        if volume < VOLUME_TOLERANCE:
            print(f"DEBUG(Volume): 计算出的体积非常小或为零 ({volume:.2e} Å^3)，请检查晶胞数据: '{lattice_str}'")
            return 0.0 # 返回0表示实际体积为零
        
        return volume
    except (ValueError, np.linalg.LinAlgError) as e:
        print(f"ERROR(Volume): 计算晶胞体积出错，晶胞字符串: '{lattice_str}'。错误: {e}")
        return 0.0 # 计算错误时返回0

def extract_stress_and_volume_from_extxyz(extxyz_file_path):
    """
    读取 .extxyz 文件，提取每帧的应力张量（9分量，eV/Å^3）和晶胞体积（Å^3）。

    Args:
        extxyz_file_path (str): .extxyz 文件的路径。

    Returns:
        tuple: (list_of_stresses, list_of_volumes)
               list_of_stresses: 列表，每个元素是该帧的 9 个应力分量组成的 NumPy 数组。
               list_of_volumes: 列表，每个元素是该帧的晶胞体积。
               如果文件无法读取或解析，返回 ([], [])。
    """
    all_stresses = []
    all_volumes = []

    if not os.path.exists(extxyz_file_path):
        print(f"ERROR: .extxyz 文件未找到: {extxyz_file_path}")
        return [], []

    try:
        # 使用 ase.io.read 来健壮地处理 extxyz 文件
        ase_atoms_list = ase.io.read(extxyz_file_path, index=':')
        print(f"DEBUG: 使用 ase.io.read 成功读取到 {len(ase_atoms_list)} 帧。")
        
        if not ase_atoms_list:
            print(f"WARNING: 在 {extxyz_file_path} 中未找到任何帧。")
            return [], []

        for i, atoms in enumerate(ase_atoms_list):
            current_stress = np.zeros(9) # 默认值，如果未找到应力则为0
            current_volume = 0.0          # 默认值，如果未找到晶胞或计算失败则为0

            # --- 提取 stress 信息 ---
            if atoms.has('stress'):
                # 获取 3x3 矩阵形式的应力，然后展平
                stress_matrix = atoms.get_stress(voigt=False) # 获取 3x3 矩阵
                
                # 如果是 Voigt 形式 (6个分量)，则转换为 3x3 矩阵
                if stress_matrix.shape == (6,):
                    sxx, syy, sizz, sxy, sxz, syz = stress_matrix
                    stress_matrix = np.array([
                        [sxx, sxy, sxz],
                        [sxy, syy, syz],
                        [sxz, syz, sizz]
                    ])
                
                # 假设单位已经是 eV/A3 (根据你的cp2k2xyz.py脚本输出)
                current_stress = stress_matrix.flatten() # 展平为 9 个分量
                print(f"DEBUG(Frame {i+1}): 成功提取应力。")
            else:
                print(f"WARNING(Frame {i+1}): 在 {extxyz_file_path} 中的当前帧未找到 'stress' 属性。将使用零填充。")

            # --- 计算晶胞体积 ---
            if atoms.has('cell') and atoms.get_cell().any():
                # atoms.get_cell().array 返回 3x3 晶胞矩阵，扁平化后传递给函数
                lattice_str = " ".join(map(str, atoms.get_cell().array.flatten()))
                temp_volume = calculate_cell_volume(lattice_str)
                if temp_volume is not None: # calculate_cell_volume 现在只返回 float (0.0) 或 int (0)
                    current_volume = temp_volume
                else: # 如果 calculate_cell_volume 内部出现严重错误导致返回 None (理论上不应该，因为我已改为返回0.0)
                    print(f"ERROR(Frame {i+1}): 晶胞体积计算函数返回 None。将使用零体积。")
            else:
                print(f"WARNING(Frame {i+1}): 在 {extxyz_file_path} 中的当前帧未找到 'cell' (晶胞) 属性或其无效。将使用零体积。")

            all_stresses.append(current_stress)
            all_volumes.append(current_volume)
            print(f"DEBUG(Frame {i+1}): stress 列表当前长度: {len(all_stresses)}, volume 列表当前长度: {len(all_volumes)}")

    except Exception as e:
        print(f"FATAL ERROR: 处理文件 {extxyz_file_path} 时发生异常: {e}")
        # 如果是这里捕获的异常，那么 all_stresses 和 all_volumes 可能是不完整的
        return [], []

    # 最终检查列表长度是否一致
    if len(all_stresses) != len(all_volumes):
        print(f"ERROR: 应力数据 ({len(all_stresses)} 帧) 和体积数据 ({len(all_volumes)} 帧) 的帧数不匹配。这不应该发生，因为每个循环都添加了数据。")
        print("请检查你的 extxyz 文件是否有异常结构或残缺帧。")
        return [], []
    
    print(f"DEBUG: 成功提取所有帧的数据。总帧数: {len(all_stresses)}")
    return all_stresses, all_volumes

def stress_datafromextxyz():
    """
    扫描当前目录下所有*.extxyz文件，提取每帧的stress（9分量，eV/A^3）和Lattice（计算体积），返回两个列表（每帧一个元素）。
    适配头行同时包含Lattice和stress的extxyz格式。
    返回: (all_stresses, all_volumes)
    """
    all_stresses = []
    all_volumes = []
    extxyz_files = glob.glob("*.extxyz")
    for extxyz_file in extxyz_files:
        with open(extxyz_file, 'r') as f:
            lines = f.readlines()
        for line in lines:
            if 'Lattice=' in line and 'stress=' in line:
                # 提取Lattice
                lattice_match = re.search(r'Lattice="([\d\.\-eE\s]+)"', line)
                if lattice_match:
                    lattice_str = lattice_match.group(1)
                    coords = list(map(float, lattice_str.split()))
                    if len(coords) == 9:
                        cell_matrix = np.array(coords).reshape(3, 3)
                        volume = abs(np.linalg.det(cell_matrix))
                    else:
                        volume = 0.0
                else:
                    volume = 0.0
                # 提取stress
                stress_match = re.search(r'stress="([\d\.\-eE\s]+)"', line)
                if stress_match:
                    stress_str = stress_match.group(1)
                    stress = np.array(list(map(float, stress_str.split())))
                    if stress.shape[0] == 9:
                        all_stresses.append(stress)
                        all_volumes.append(volume)
                    else:
                        all_stresses.append(np.zeros(9))
                        all_volumes.append(volume)
                else:
                    all_stresses.append(np.zeros(9))
                    all_volumes.append(volume)
    return all_stresses, all_volumes

def calculate_virial_from_stress(stress_tensor, volume):
    """
    根据应力张量、体积和给定公式计算维里张量。
    """
    if volume == 0.0: # Changed from 'is None or volume == 0.0'
        # print(f"Warning: Cannot calculate virial with zero volume. Stress: {stress_tensor}")
        return np.zeros_like(stress_tensor) # Return zeros if volume is problematic
        
    if not isinstance(stress_tensor, np.ndarray) or stress_tensor.shape != (9,):
        print(f"ERROR: 提供的应力张量格式无效: {stress_tensor}")
        return None

    # 您提供的转换公式：stress * 体积 * 10^4 / (6.24 * 1602.17)
    conversion_factor = -10**4 / (6.24 * 1602.17)
    
    virial_tensor = stress_tensor * volume * conversion_factor
    
    return virial_tensor

def write_virial_raw(output_file_path, all_virials):
    """
    将维里张量数据写入 virial.raw 文件。
    文件中的数据是扁平化的。
    """
    if not all_virials:
        print("WARNING: 没有维里数据可写入。跳过文件创建。")
        return

    flat_virials = np.concatenate([v.flatten() for v in all_virials])

    try:
        np.savetxt(output_file_path, flat_virials, fmt='%20.10f')
        print(f"成功将维里数据写入: {output_file_path}")
    except Exception as e:
        print(f"ERROR: 写入维里数据到 {output_file_path} 时出错: {e}")

def write_virial_raw_from_stress_data(stresses, volumes, output_file="virial.raw",outdir="./data"):
    """
    根据stresses和volumes，按公式输出九列数据到virial.raw，首行为空。
    """
    if not stresses or not volumes or len(stresses) != len(volumes):
        print("输入数据为空或长度不一致，无法写入virial.raw")
        return
    conversion_factor = 1e4 / (6.24 * 1602.17)
    os.makedirs(outdir, exist_ok=True)
    output_path = os.path.join(outdir, output_file)
    with open(output_path, 'w') as f:
        for stress, volume in zip(stresses, volumes):
            virial = stress * volume * conversion_factor
            f.write(" ".join(f"{v:20.10f}" for v in virial) + "\n")
    print(f"已写入 {len(stresses)} 行到 {output_file}")

def write_stress_and_volume_data(stresses, volumes, outdir="."):
    """
    将stress和volume分别输出到data/stress.txt和data/volume.txt。
    每行对应一帧，stress为九列，volume为一列。
    """
    os.makedirs(outdir, exist_ok=True)
    stress_path = os.path.join(outdir, "stress.txt")
    volume_path = os.path.join(outdir, "volume.txt")
    with open(stress_path, 'w') as f_stress:
        for stress in stresses:
            f_stress.write(" ".join(f"{v:20.10f}" for v in stress) + "\n")
    with open(volume_path, 'w') as f_volume:
        for v in volumes:
            f_volume.write(f"{v:20.10f}\n")
    print(f"已写入 {len(stresses)} 行到 {stress_path} 和 {volume_path}")

# --- 主执行部分 (保持不变，用于测试) ---
if __name__ == "__main__":
    # 1. 模拟一个 extxyz 文件用于测试
    # 注意：这个模拟文件需要包含 stress 和 Lattice 信息
    mock_extxyz_content = """
3
energy=-123.456 pbc="T T T" Lattice="10.0 0.0 0.0 0.0 10.0 0.0 0.0 0.0 10.0" Properties=species:S:1:pos:R:3:stress:R:9
H    0.0 0.0 0.0   stress="1.0 0.1 0.2 0.3 2.0 0.4 0.5 0.6 3.0"
O    1.0 0.0 0.0
H    0.0 1.0 0.0
3
energy=-124.567 pbc="T T T" Lattice="10.1 0.0 0.0 0.0 10.1 0.0 0.0 0.0 10.1" Properties=species:S:1:pos:R:3:stress:R:9
H    0.1 0.0 0.0   stress="1.1 0.11 0.22 0.33 2.1 0.44 0.55 0.66 3.1"
O    1.1 0.0 0.0
H    0.1 1.0 0.0
"""
    test_extxyz_file = "test_stress.extxyz"
    with open(test_extxyz_file, "w") as f:
        f.write(mock_extxyz_content)
    print(f"Generated mock .extxyz file: {test_extxyz_file}\n")

    # 2. 提取应力张量和晶胞体积
    print("--- Extracting stress and volume from .extxyz ---")
    stresses_from_extxyz, volumes_from_extxyz = extract_stress_and_volume_from_extxyz(test_extxyz_file)

    if stresses_from_extxyz and volumes_from_extxyz:
        print(f"\n成功提取到 {len(stresses_from_extxyz)} 帧应力和体积数据。")
        print(f"第一帧应力 (eV/Å^3):\n{stresses_from_extxyz[0]}")
        print(f"第一帧体积 (Å^3): {volumes_from_extxyz[0]:.6f}\n")

        # 3. 计算所有帧的维里张量 (计算公式和写入部分不变)
        print("--- Calculating virial tensors ---")
        all_calculated_virials = []
        for i in range(len(stresses_from_extxyz)):
            virial = calculate_virial_from_stress(stresses_from_extxyz[i], volumes_from_extxyz[i])
            if virial is not None:
                all_calculated_virials.append(virial)
        
        if all_calculated_virials:
            print(f"成功计算到 {len(all_calculated_virials)} 帧维里数据。")
            print(f"第一帧计算出的维里 (eV):\n{all_calculated_virials[0]}\n")

            # 4. 将维里数据写入 virial.raw 文件
            output_raw_file = "virial.raw"
            write_virial_raw(output_raw_file, all_calculated_virials)
        else:
            print("未计算任何维里数据。")
    else:
        print("未能提取应力和体积数据。")

    # 新增批量处理extxyz文件的功能
    stresses, volumes = stress_datafromextxyz()
    write_virial_raw_from_stress_data(stresses, volumes)
    write_stress_and_volume_data(stresses, volumes)

    # 清理测试文件 (根据需要启用或禁用)
    # os.remove(test_extxyz_file)
    # if os.path.exists("virial.raw"):
    #     os.remove("virial.raw")

