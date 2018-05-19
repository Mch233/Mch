from sklearn.datasets import make_blobs
from matplotlib import pyplot
import numpy as np


#########################################数据生成################################################
n = input('输入样本个数\n') 

#n = int(n)-3
#center=np.array([[-1,1],[1,1],[2,2]])#中心点
#labelcen=np.array([0.0,1.0,2.0])
n = int(n)-2
center=np.array([[3,8],[8,4]])#中心点
labelcen=np.array([0.0,1.0])

datacen = np.c_[center,labelcen]

data, label = make_blobs(n_samples=n, n_features=2, centers=center)
alldata = np.c_[data,label]#合并xy和label

alldata = np.r_[alldata,datacen]
print(center)

#########################################绘制样本显示############################################


pyplot.scatter(data[:, 0], data[:, 1], c=label)

pyplot.scatter(center[:, 0], center[:, 1], marker="*",c=[0.0,1.0])#中心点

pyplot.show()

###########################################写入文件#############################################

np.savetxt('data.txt',alldata,fmt='%.2f')
        
#np.savetxt('data.txt',alldata,fmt='%.2f',delimiter=' ',newline='\n' )
