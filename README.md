# scripts
# CMS Scripts
计算材料科学相关代码/脚本  
**注意事项：**  
• 本仓库代码主要依赖dpdata、ase、pymatgen等python包    
• 使用脚本前，请先阅读源码及其中的注释！  
• 部分脚本可解析命令行参数，使用 `python xx.py -h`  

---
# 脚本内容  




Markdown 书写指令基本需求  

[Gemini链接](https://gemini.google.com/app/25f8fe3324929ac5)）  
![Markdown Logo](https://markdown.com.cn/images/logo.png)  
` print(“基本用法”）`

## VASP相关脚本
• VASP-scripts/:DeepMD相关脚本
```python
outcar_vasp2raw.py                   #将OUTCAR转换为DeePMD的输入文件（利用 dpdata package）
outcar_vasp2extxyz.py                #将OUTCAR转换为MACE的输入文件 (利用 ase package)

```



## plot相关脚本
• CP2K-script/:画图相关脚本
```python
cp2k2extxyz.py                       #将cp2k-aimd文件转换为extxyz文件（基于ase以及brucefan1983GPUMD)
vasp_aimd_plot.sh                    #提取VASP的AIMD输出文件信息绘图 (基于matplotlib) 包含能量、体积、温度、应力等
```


## lammps相关脚本
• Lammps-scripts/Lammps相关脚本
```python
run-msd-mace.in                      #运行lammps的NVT系统MD，输出msd数据存储到msd_results.txt;
```



## plot相关脚本
• plots/:画图相关脚本
```python
plot-dp-test.py                      #直接对deepmd测试的输出文件绘图（基于matplotlib)
vasp_aimd_plot.sh                    #提取VASP的AIMD输出文件信息绘图 (基于matplotlib) 包含能量、体积、温度、应力等
```



## strucs相关脚本
• strucs-scripts/结构处理:的相关脚本
```python
check_density.py                     #处理结构的密度（利用 ase package python check_axis.py  *.vasp)
check_axis.py                        #处理结构的abc以及aplha beta gamma晶格参数 执行方式python check_axis.py  *.vasp)
check_nd.py                          #处理结构得到中子衍射数据(利用 deby package) 
change_poscar2xyz.py                 #转换结构形式，获得xyz格式的文件
```

## 性质计算脚本
• properties-scripts/性质计算的相关脚本
```python
calculate_micro_ion_conductivity.py  #微观角度计算离子电导率的计算脚本（ 利用ase package)
```



## Machinelearning相关脚本
• machinelearning-scripts/机器学习的相关脚本
```python
chgnet_train.py			    #chgnet训练脚本
chgnet_grep.py			    #chgnet提取vasp的输出文件
```
