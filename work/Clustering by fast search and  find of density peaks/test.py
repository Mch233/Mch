from numpy import *
import numpy as np
import matplotlib.pyplot as plt
import math

file0 = input('数据\n')
f = loadtxt(file0)
plt.scatter(f[:,0],f[:,1],c = f[:,2])
plt.show()
      

