#!/bin/bash
# example file is cp2k_aimd_2025.1 
cp *.out output.log
python 01_cp2k2extxyz.py
python 02_cp2k_2npy_except_stress.py
python 03_extxyz2stress_raw.py 
python 04_cp2k_2npy.py
cp original-stress.extxyz mace_trainsey.extxyz
