import sys
from ase.io import read
from ase import units

def calculate_density(atoms):
    cell_volume = atoms.get_volume() * 1e-24 # 从 Å^3 转换为 cm^3
    mass = atoms.get_masses().sum() * units._amu # 将原子质量单位转换为kg（ASE中的amu是以kg为单位的）
    density = mass / cell_volume # 密度单位是 kg/cm^3
    density *= 1e3  # 将密度从 kg/cm^3 转换为 g/cm^3
    return density

if __name__ == "__main__":
    if len(sys.argv) != 2:
        print("Usage: python density.py POSCAR")
        sys.exit(1)

    poscar_path = sys.argv[1]

    # 读取POSCAR文件
    atoms = read(poscar_path)

    # 计算晶胞密度
    density = calculate_density(atoms)

    # 打印结果
    print(f"The density of the crystal is {density:.3f} g/cm^3")

