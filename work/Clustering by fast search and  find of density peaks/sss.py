from matplotlib import pyplot
import numpy as np
import matplotlib.pyplot as plt    
from sklearn.datasets import make_gaussian_quantiles  

#########################################数据生成################################################

X1, Y1 = make_gaussian_quantiles(n_samples=1000,n_features=2, n_classes=4)  

plt.scatter(X1[:, 0], X1[:, 1], marker='o', c=Y1) 
alldata = np.c_[X1,Y1]



###########################################写入文件#############################################

np.savetxt('sss.txt',alldata,fmt='%.2f')

#np.savetxt('data.txt',alldata,fmt='%.2f',delimiter=' ',newline='\n' )
