from numpy import *
import numpy as np
import matplotlib.pyplot as plt
import math
from sklearn.datasets import make_gaussian_quantiles 

#########################################打开文件###############################################
filein = input('输入原始数据文件\n')
f = loadtxt(filein)
N = len(f)

############################################文件转换#############################################
fileout = input('输入实验数据文件\n')
b = open(fileout,'w')

for i in range(1,N):
    for j in range(i+1,N+1):
        x = float(f[i-1][0]) - float(f[j-1][0])
        y = float(f[i-1][1]) - float(f[j-1][1])
        length = math.sqrt(x*x + y*y)
        length = ("%.2f" %length)
        
        b.write(str(i))
        b.write(' ')
        b.write(str(j))
        b.write(' ')
        b.write(str(length))
        b.write('\n')

b.close()

