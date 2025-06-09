# 这个脚本用于最终计算离子电导率，微观角度进行计算，最终数值基本一致与宏观载流子电导率计算  
# 计算离子电导率
# 需要的参数：count, slope, dim_fac, cell_volume, temperature
# count: Li离子数目
# slope: 从MSD计算得到的斜率
# dim_fac: 维度因子，通常为3
# cell_volume: 单位胞体积，单位为A^3
# temperature: 温度，单位为K
# 计算离子电导率
# 公式：dd = slope / (2 * dim_fac) * 1e-4 * conversion_cm
# conversion_cm: 转换因子，将结果从A^2/ps转换为mS/cm
# 公式：dd = slope / (2 * dim_fac) * 1e-4 * conversion_cm
# conversion_cm = 1000 * count / (cell_volume * 1e-24 * const.N_A) * (const.N_A * const.e)**2 / (const.R * temperature)
#得到dd= mS/cm
import ase.constants as const
from scipy.stats import linregress
import ase.traj
# 计算离子电导率
    conversion_cm = 1000 * count / \
        (cell_volume * 1e-24 * const.N_A) * \
        (const.N_A * const.e)**2 / (const.R * temperature)
dd=slope / (2*MSD.dim_fac) * 1e-4*conversion_cm

#ount=96, slope这里是0.02，dim_fac=3，体积是8852A3，温度400K，计算dd,dd=16.8mS/cm
