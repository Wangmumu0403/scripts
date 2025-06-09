#!/bin/bash

# --- 1. 数据提取 (Bash 命令) ---
# 这些命令将从 VASP 输出文件中提取所需的数据，并保存到独立的文件中。
# 使用 '>>' (追加) 操作符意味着如果文件已存在，新数据将添加到文件末尾。
# 在运行此脚本之前，请确保 'OUTCAR' 和 'OSZICAR' 文件存在于当前目录。
#
# 注意：
# 你提供的 grep 命令 (例如 'grep "F=" OSZICAR >> energy') 将直接把匹配到的行
# 追加到指定文件。你的 Python 脚本期望 'energy' 文件有特定列（0, 2, 4），
# 'volume' 文件有第 4 列，'stress' 文件有 2-7 列。
# 请务必手动检查这些 grep 命令生成的文件的实际列数和内容是否与你的 Python 脚本
# 中 np.loadtxt 的 'usecols' 参数精确匹配。如果它们不匹配，Python 脚本可能会报错。

#echo "正在从 OUTCAR 和 OSZICAR 文件中提取数据..."

# 提取能量数据 (从 OSZICAR 的 "F=" 行)
grep "F=" OSZICAR >> energy

# 提取体积数据 (从 OUTCAR 的 "volume" 行)
grep "volume" OUTCAR >> volume

# 提取应力数据 (从 OUTCAR 的 "in kB" 行)
grep "in kB" OUTCAR >> stress

echo "数据提取完成。正在运行 Python 绘图脚本..."

# --- 2. Python 绘图脚本 ---
# 以下是你的 Python 脚本，它将读取上述生成的数据文件并绘制图表。
# 请将此整个 Python 代码块保存为一个独立的 .py 文件 (例如：'plot_vasp_npt.py')
# 然后在下面使用 'python3 plot_vasp_npt.py' 来执行。

python3 << EOF
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker # Import the ticker module

# --- 1. Data Loading ---
# Load your VASP data files using the specified np.loadtxt commands

try:
    # Energy data: Assumes 'energy' file has ionic step in column 0,
    # TEMPERATURE in column 2 (index 2), and total energy (F=) in column 4 (index 4).
    # Please verify these column indices with your actual 'energy' file content!
    data0 = np.loadtxt("energy", usecols=(0, 2, 4), comments="#") # Ionic Step, Temperature, Total Energy (F=)
except FileNotFoundError:
    print("Error: 'energy' file not found. Ensure it exists and is correctly generated.")
    exit()
except Exception as e:
    print(f"Error reading 'energy' file: {e}")
    exit()

try:
    # Volume data: Assumes 'volume' file has volume in column 4.
    data1 = np.loadtxt("volume", usecols=(4), comments="#")
except FileNotFoundError:
    print("Error: 'volume' file not found. Ensure it exists and is correctly generated.")
    exit()
except Exception as e:
    print(f"Error reading 'volume' file: {e}")
    exit()

try:
    # Stress data: Loading all six stress components (columns 2, 3, 4, 5, 6, 7).
    data2 = np.loadtxt("stress", usecols=(2, 3, 4, 5, 6, 7), comments="#")
except FileNotFoundError:
    print("Error: 'stress' file not found. Ensure it exists and is correctly generated.")
    exit()
except Exception as e:
    print(f"Error reading 'stress' file: {e}")
    exit()

# --- 2. Create a Single Figure with Subplots ---

fig, axes = plt.subplots(4, 1, figsize=(12, 22)) # 4 rows, 1 column of plots for more space

# Plot 1: Energy Convergence (axes[0])
axes[0].plot(data0[:, 0], data0[:, 2], 'o-', markersize=4, label='Total Energy (F=)') # data0[:, 2] is F= Energy
axes[0].set_xlabel('Ionic Step')
axes[0].set_ylabel('Total Energy (eV)')
axes[0].set_title('Energy Convergence')
axes[0].legend()
axes[0].grid(True)
axes[0].yaxis.set_major_formatter(ticker.FormatStrFormatter('%.4f')) # Formats to 4 decimal places

# Plot 2: All Six Stress Components Evolution (axes[1])
stress_labels = ['Stress_XX', 'Stress_YY', 'Stress_ZZ', 'Stress_XY', 'Stress_YZ', 'Stress_ZX']
for i in range(data2.shape[1]): # Iterate through all 6 columns of stress data
    axes[1].plot(np.arange(data2.shape[0]), data2[:, i], label=stress_labels[i])
axes[1].set_xlabel('Simulation Step')
axes[1].set_ylabel('Stress (kB)')
axes[1].set_title('Stress Components Evolution')
axes[1].legend(ncol=2) # Arrange legend in two columns for better fit
axes[1].grid(True)

# Plot 3: Temperature Evolution vs. Ionic Step (NEW POSITION: axes[2])
# data0[:, 0] is Ionic Step, data0[:, 1] is Temperature
axes[2].plot(data0[:, 0], data0[:, 1], 'o-', markersize=4, label='Temperature')
axes[2].set_xlabel('Ionic Step')
axes[2].set_ylabel('Temperature (K)')
axes[2].set_title('Temperature Evolution')
axes[2].legend()
axes[2].grid(True)

# Plot 4: Volume Evolution vs. Ionic Step (NEW POSITION: axes[3])
min_len_vol = min(data0.shape[0], len(data1))
axes[3].plot(data0[:min_len_vol, 0], data1[:min_len_vol], 'o-', markersize=4, label='Volume')
axes[3].set_xlabel('Ionic Step')
axes[3].set_ylabel('Volume (Å$^3$)')
axes[3].set_title('Cell Volume Evolution')
axes[3].legend()
axes[3].grid(True)


# Adjust layout to prevent labels from overlapping
plt.tight_layout()

# Display the combined plot
plt.show()
plt.savefig("AIMD-NPT-energy.png", dpi=300, bbox_inches="tight")

EOF

echo "绘图完成。结果已保存到 AIMD-NPT-energy.png"
