from numpy import *
import numpy as np
import matplotlib.pyplot as plt

#读取文件
test = input('输入\'文件名\'\n')
f = loadtxt(test)

#计算数据点个数
N = max(f[:,0])
D = max(f[:,1])
ND = N
if N < D:
    ND = D
ND = int(ND)

#向矩阵中装载数据
A = zeros((ND,ND),dtype=float)
for line in f:
    row=int(line[0]) - 1
    col=int(line[1]) - 1
    A [row,col] = float(line[2])
    A [col,row] = float(line[2])
print('矩阵A')
print(A)

#计算bc
val=f[:,2]
N = len(val)
t = 2.0
site = round(t/100*N)
val = sort(val)
bc = val[site -1 ]
print('bc =')
print(bc)

#计算密度rho
#rho = [0]*ND
rho = np.zeros(ND)
for raw in range(ND):
    for col in range(ND):
        if raw != col:
            if A[raw,col] < bc:
                rho[raw] = rho[raw] + 1
print('密度rho')
print(rho)

#计算delta
key = [-1]*ND
delta = np.zeros(ND)
index =np.argsort(-rho)#按rho降排，返回索引的list
for i in index:
    ii = 1#标记ｉ的索引
    if i != index[0]:
        temp = A[i,:]#剩余局部最大密度所在行  
        mindelta = 99999.0
        j = 0
        while j < ii:#该区间索引对应index存放的密度更大点
            if temp[index[j]] < mindelta:
                mindelta = temp[index[j]]
                key[i] = index[j]
            j = j + 1
        delta[i] = mindelta
    ii = ii + 1#追踪ｉ索引的变化
delta[index[0]] = max(A[index[0],:])
print('delta')
print(delta)
print('delta对应点')
print(key)

#决策图
plt.title('decision diagram')
plt.xlabel(r'$\rho$')
plt.ylabel(r'$\delta$')
x = rho
y = delta
plt.plot(x,y,'o')
plt.show()


#gama图
gama = np.zeros(ND)
gama = rho * delta
print('gama')
print(gama)
index = np.argsort(-gama)
gama = sorted(gama,reverse=True)
print('降序gama')
print(gama)

plt.title('gama')
plt.xlabel('n')
plt.ylabel('gama')
x = range(ND)
y = gama
plt.plot(x,y,'o')
plt.show()

#选择中心个数n
n = input('输入中心点个数\n')
n = int(n)

#划分中心点
label = [-1]*ND
print('未聚类')
print(label)
k = 1
i = 0
while k <= n:
    label[index[i]] = k
    k = k + 1
    i = i + 1

print('确定中心')
print(label)


#聚类
i = 0
while i < ND:
    if label[i] == -1:
        label[i] =label[key[i]] 
    i = i + 1

print('初步聚类')
print(label)

#噪点
bord_rho = [0]*n#存放截止密度，索引=类别标号-1
for i in range(0,ND-1):
    for j in range(i+1,ND):
        
        if A[i][j] < bc and label[i] != label[j]:
            rho_aver = (rho[i]+rho[j])/2
        
            if rho_aver > bord_rho[label[i]-1]:
                bord_rho[label[i]-1] = rho_aver
            if rho_aver > bord_rho[label[j]-1]:
                bord_rho[label[j]-1] = rho_aver

print('截止密度')
print(bord_rho)

for i in range(ND):
    if rho[i] <  bord_rho[label[i]-1]:
        label[i] = 0

print('聚类结果')
print(label)


