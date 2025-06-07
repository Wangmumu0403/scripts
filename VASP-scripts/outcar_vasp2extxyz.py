from ase.io import read, write
from ase.calculators.vasp import Vasp

#支持多帧结构
# 读取 OUTCAR 中的所有离子步（返回 Atoms 对象列表）
mace_trainsets = read('OUTCAR', index=':', format='vasp-out')
# 保存为多帧 extxyz 文件（支持轨迹）
write('mace_multisets.extxyz', mace_trainsets, format='extxyz')



# 支持单帧结构
# 读取 OUTCAR 文件（自动识别格式）
atoms = read('OUTCAR')
# 保存为 extxyz 格式（包含能量、力等属性）
write('mace_singleset.extxyz', atoms, format='extxyz')
