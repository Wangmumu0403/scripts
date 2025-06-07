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
outcar_vasp2raw.py          #将OUTCAR转换为DeePMD的输入文件（基于dpdata）
outcar_vasp2extxyz.py       #将OUTCAR转换为MACE的输入文件 (基于ase)

```
## Plot相关脚本
• Plot/:画图相关脚本
```python
plot-dp-test.py             #将deepmd输出文件进行画图（基于matplotlib)
```
