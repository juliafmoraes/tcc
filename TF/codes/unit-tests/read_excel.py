#!/usr/bin/env python
import numpy as np
import pandas as pd
import math
import openpyxl
import xlrd

N = pd.read_excel(r'C:\Users\Julia\Documents\tcc\TF\codes\results\20191012\slices\trap_detrap_1.xlsx', sheet_name='N')
N_dif = pd.read_excel(r'C:\Users\Julia\Documents\tcc\TF\codes\results\20191012\slices\trap_detrap_1.xlsx', sheet_name='N_dif')
N_trap = pd.read_excel(r'C:\Users\Julia\Documents\tcc\TF\codes\results\20191012\slices\trap_detrap_1.xlsx', sheet_name='N_trap')

N = N.values.tolist()