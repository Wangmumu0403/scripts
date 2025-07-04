import dpdata
import numpy as np
# 使用LabeledSystem读取CP2K的AIMD输出文件；注意第一个参数是aimd输出所在的文件夹位置
#data=dpdata.LabeledSystem('./',cp2k_output_name='output.log',fmt='cp2kdata/md')
data=dpdata.LabeledSystem('./',cp2k_output_name='output.log', ensemble_type="NVT", fmt="cp2kdata/md")

print(data)
# 转化为deepmd的数据格式并输出到指定位置
data.to_deepmd_npy('./data',fmt='deepmd-npy')
data.to_deepmd_raw('./data',fmt='deepmd-raw')
