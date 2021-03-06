#RBM

##能量函数

为什么要弄这个能量模型呢？

	第一、RBM网络是一种无监督学习的方法，无监督学习的目的是最大可能的拟合输入数据，所以学习RBM网络的目的是让RBM网络最大可能地拟合输入数据。

	第二、对于一组输入数据来说，现在还不知道它符合那个分布，那是非常难学的。例如，知道它符合高斯分布，那就可以写出似然函数，然后求解，就能求出这个是一个什么样个高斯分布；但是要是不知道它符合一个什么分布，那可是连似然函数都没法写的，问题都没有，根本就无从下手。

	任何概率分布都可以转变成基于能量的模型

	第三、在马尔科夫随机场（MRF）中能量模型主要扮演着两个作用：一、全局解的度量（目标函数）；二、能量最小时的解（各种变量对应的配置）为目标解。也就是能量模型能为无监督学习方法提供两个东西：a）目标函数；b）目标解。

对于一组给定的状态(v,h)，定义如下的能量函数：


每个可视节点和隐藏节点之间的连接结构都有一个能量


求解输入样本的极大似然，就能让RBM网络表示的Gibbs分布和样本本身表示的分布最接近

##概率


1. 联合概率
对能量函数指数化和正则化后可以得到可见层节点集合和隐藏层节点集合分别处于某一种状态下

归一化因子表示对可见层和隐藏层节点集合的所有可能状态的（能量的指数形式）求和。

2. 边缘分布（似然函数）
可见层节点集合处于某一种状态分布下的概率

通过对隐藏层节点集合的所有状态求和，可以得到可见层节点集合的边缘分布
Gibbs分布的概率密度函数
求解的目标——让RBM网络的表示Gibbs分布最大可能的拟合输入数据。

定义样本表示的分布和RBM网络表示的边缘分布的KL距离


##模型训练

v 显性神经元  （接受输入）
h 隐性神经元  （提取特征）特征探测器
连接方式：对称 (双向) 全连接

显层和隐层内部的神经元都没有互连，只有层间的神经元有对称的连接线。
这样的好处是，在给定所有显元的值的情况下，每一个隐元取什么值是互不相关的。
可并行地计算整层神经元。

每一个节点（无论是Hidden Unit还是Visible Unit）都有两种状态：
处于激活状态时值为1，未被激活状态值为0
节点的激活概率由可见层和隐藏层节点的分布函数计算。

RBM 模型需要确定两部分。
1. 如果想确定这个模型，首先是要知道可见层和隐藏层节点个数，
可见层节点个数即为输入的数据维数，
隐藏层节点个数在一些研究领域中是和可见层节点个数有关的，
多数情况下，隐藏层节点个数需要根据使用而定
或者是在参数一定的情况下，使得模型能量最小时的隐藏层节点个数。

2. 要想确定这个模型还得要知道模型的三个参数θ=W,a,b
		可见层和隐藏层之间的权重，
		可见层的偏置
		隐藏层的偏置

对于给定的训练样本，通过训练得到参数θ，使得在该参数下，由RBM表示的概率分布尽可能与训练数据相符合。
最大化似然函数;			常用的方法是梯度下降法

RBM应用到了能量学上的一些知识：在能量最少的时候，物质最稳定。
应用到RBM就是，在能量最少的时候，网络最稳定，也就是网络最优。
能量E和概率P是成反比的关系，所以通过最大化P，才能使能量值E最小
因此此时的损失函数可以看做是-P(x)，且求导时需要是加上负号的。


##RBM用途
RBM的用途主要是两种，
1. 对数据进行编码，然后交给监督学习方法去进行分类或回归，
2. 得到了权重矩阵和偏移量，供BP神经网络初始化训练。
