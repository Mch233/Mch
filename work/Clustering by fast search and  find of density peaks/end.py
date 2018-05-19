from numpy import *
import numpy as np
import matplotlib.pyplot as plt
import math

##############################################读取文件#############################################
test = input('输入\'文件名\'\n')
f = loadtxt(test)

#########################################计算数据点个数############################################
N = max(f[:,0])
D = max(f[:,1])
ND = N
if D > N:
    ND = D
ND = int(ND)
print (ND)

##########################################向矩阵中装载数据#########################################
A = np.zeros((ND,ND),dtype=float)
for line in f:
    row=int(line[0]) - 1
    col=int(line[1]) - 1
    A [row,col] = float(line[2])
    A [col,row] = float(line[2])
print('矩阵A')
print(A)
np.savetxt('Matrix',A,fmt='%.2f',newline='\n',delimiter=' ')

#############################################计算bc############################################
val=f[:,2]
N = len(val)
t = 2.0
site = round(t/100*N)
val = sort(val)
bc = val[site -1 ]
print('bc =')
print(bc)

##############################################计算密度rho######################################
rho = np.array([0.0]*ND)
for raw in range(ND-1):
    for col in range(raw+1,ND):
                t = math.exp(-(A[raw,col]/bc)*(A[raw,col]/bc))
                rho[raw] = rho[raw] + t
                rho[col] = rho[col] + t
print('密度rho')
print(rho)

###############################################计算delta#####################################
key =np.array([-1.0]*ND)
delta = np.array([0.0]*ND)
index =np.argsort(-rho)#按rho降排，返回索引的list
for i in range(1,len(index)):
    temp = A[index[i],:]
    mindelta = 99999.9
    j = 0
    while j < i:
        if mindelta > temp[index[j]]:
            mindelta = temp[index[j]]
            delta[index[i]] = mindelta
            key[index[i]] = index[j]
        j = j + 1
delta[index[0]] = max(A[index[0],:])

print('delta')
print(delta)
print('delta对应点')
print(key)

########################################写入delta rho#######################################
filerd = open('filerd','w')
for i in range(ND):
    filerd.write(str(rho[i]))
    filerd.write(' ')
    filerd.write(str(delta[i]))
    filerd.write('\n')

##############################################决策图#########################################
plt.title('decision diagram')
plt.xlabel(r'$\rho$')
plt.ylabel(r'$\delta$')
x = rho
y = delta
plt.plot(x,y,'o')
plt.show()


################################################gama图#####################################
gama = np.zeros(ND)
gama = rho * delta
print('gama')
print(gama)
ind = np.argsort(-gama)
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

###############################################选择中心个数n################################
n = input('输入中心点个数\n')
n = int(n)

#################################################划分中心点################################
label = [-1]*ND
central = [0]*n
print('未聚类')
print(label)
k = 1
i = 0
while k <= n:
    if delta[ind[i]] > bc:
        label[ind[i]] = k
        central[k-1]=ind[i]
        k = k + 1
    i = i + 1
 
print('确定中心')
print(label)
print(central)


##############################################聚类###########################################
i = 0
while i < ND:
    if label[index[i]] == -1:
        label[index[i]] = label[int(key[index[i]])] 
    i = i + 1

label0 = label
print('初步聚类')
print(label0)

##############################################噪点###########################################
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
    if rho[i] <  bord_rho[label[i]-1]:#该点的密度是否小于对应类标的截止密度
        #label[i]类标;labe[i]-1类标在截止密度list中的索引
        label[i] = 0
label1 = label
print('聚类结果')
print(label1)
#################################################################################################
#打开文件
filein = input('输入原始数据文件\n')
f = loadtxt(filein)
N = len(f)

#############################################数据分布图#######################################
def listofloat(L=[]):
    for i in range(len(L)):
        L[i] = float(L[i])
    return L

x = listofloat(f[:,0])
y = listofloat(f[:,1])

##############################################label0#######################################
plt.title('start')
plt.xlabel('x')
plt.ylabel('y')
colors = label0
plt.grid(True)#标尺
plt.scatter(x, y, c=colors)
plt.show()

##############################################label#######################################
plt.title('end')
plt.xlabel('x')
plt.ylabel('y')
colors = label1
plt.grid(True)#标尺
plt.scatter(x, y, c=colors)
plt.show()

