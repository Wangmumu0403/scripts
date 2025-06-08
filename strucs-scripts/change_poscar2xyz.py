import dpdata
  
# -----------------------------------
# 示例 1: 转换单个结构（如 POSCAR）
# -----------------------------------
# 读取 POSCAR
system = dpdata.System("POSCAR", fmt="vasp/poscar")

# 保存为 XYZ
system.to("xyz", "output.xyz")
