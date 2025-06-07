import dpdata
  
# 读取OUTCAR文件
system = dpdata.LabeledSystem("OUTCAR", fmt="vasp/outcar")

# 保存为DeePMD的raw格式（自动生成type.raw、coord.raw等）
system.to_deepmd_raw("deepmd_data")
system.to_deepmd_npy("deepmd_data")  # 可选：保存为npy格式
