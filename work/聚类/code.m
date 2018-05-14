clear all
close all

%% 从文件中读取数据
mdist=input('name of the distance matrix file (with single quotes)?\n');
xx=load(mdist);

ND=max(xx(:,2));　
NL=max(xx(:,1));
if (NL>ND)
  ND=NL;    %% 确保 DN 取为第一二列最大值中的较大者，并将其作为数据点总数
end

N=size(xx,1); %% xx 第一个维度的长度，相当于文件的行数（即距离的总个数）

%% 初始化矩阵dist　ND*ND
for i=1:ND
  for j=1:ND
    dist(i,j)=0;
  end
end

%% 构建距离矩阵
%% 不考虑对角线元素
for i=1:N
  ii=xx(i,1);
  jj=xx(i,2);
  dist(ii,jj)=xx(i,3); 　   ％％上三角　
  dist(jj,ii)=xx(i,3);　　％％下三角
end

%% 确定 dc
percent=2.0;
fprintf('average percentage of neighbours (hard coded): %5.6f\n', percent);
position=round(N*percent/100);  %% round 是一个四舍五入函数
sda=sort(xx(:,3));    %% 对所有距离值作升序排列
dc=sda(position);

%% 计算局部密度 rho (利用 Gaussian 核)
fprintf('Computing Rho with gaussian kernel of radius: %12.6f\n', dc);


％％计算密度
%% 数据的点密度向量初始化为０
for i=1:ND
  rho(i)=0.;
end

for i=1:ND-1　
  for j=i+1:ND
     rho(i)=rho(i)+exp(-(dist(i,j)/dc)*(dist(i,j)/dc));
     rho(j)=rho(j)+exp(-(dist(i,j)/dc)*(dist(i,j)/dc));
  end
end

%
% "Cut off" kernel
%
%for i=1:ND-1
%  for j=i+1:ND
%    if (dist(i,j)<dc)%%对称结构　ij,ji=1
%       rho(i)=rho(i)+1.;
%       rho(j)=rho(j)+1.;
%    end
%  end
%end


%%计算delta
%% 先求矩阵的列最大值，再求最大值，最后得到所有距离值中的最大值
maxd=max(max(dist));

%% 将 rho 按列降序排列，ordrho 原序
[rho_sorted,ordrho]=sort(rho,'descend');

%% 处理 rho 值最大的数据点
delta(ordrho(1))= -1.;                      ％％ordrho(1) rho最高对应的ｉ
nneigh(ordrho(1))=0;                        ％％距离

%% 生成 delta 和 nneigh 数组
for ii=2:ND　　　  ％％依次求最短距离
   delta(ordrho(ii))=maxd;
   for jj=1:ii-1　　　％％密度更大数据点的
     if(dist(ordrho(ii),ordrho(jj))<delta(ordrho(ii)))
            delta(ordrho(ii))=dist(ordrho(ii),ordrho(jj));　　％％rho 值更大的距离最近点
            nneigh(ordrho(ii))=ordrho(jj);　　％％与 ordrho(ii) 距离最近的点的编号
        %% 记录 rho 值更大的数据点中与 ordrho(ii) 距离最近的点的编号 ordrho(
     end
   end
end

%% 生成 rho 值最大数据点的 delta 值
delta(ordrho(1))=max(delta(:));


%% 决策图
disp('Generated file:DECISION GRAPH')
disp('column 1:Density')
disp('column 2:Delta')

％％保存　rho,delta
fid = fopen('DECISION_GRAPH', 'w');
for i=1:ND
   fprintf(fid, '%6.2f %6.2f\n', rho(i),delta(i));
end

%% 选择一个围住类中心的矩形
disp('Select a rectangle enclosing cluster centers')
%% 每台计算机，句柄的根对象只有一个，就是屏幕，它的句柄总是 0
%% >> scrsz = get(0,'ScreenSize')
%% scrsz =
%%            1           1        1280         800
%% 1280 和 800 就是你设置的计算机的分辨率，scrsz(4) 就是 800，scrsz(3) 就是 1280
scrsz = get(0,'ScreenSize');
%% 人为指定一个位置，感觉就没有那么 auto 了 :-)
figure('Position',[6 72 scrsz(3)/4. scrsz(4)/1.3]);

％％ gamma
for i=1:ND
  ind(i)=i;
  gamma(i)=rho(i)*delta(i);
end

%% 利用 rho 和 delta 画出一个“决策图”

subplot(2,1,1)
tt=plot(rho(:),delta(:),'o','MarkerSize',5,'MarkerFaceColor','k','MarkerEdgeColor','k');
title ('Decision Graph','FontSize',15.0)
xlabel ('\rho')
ylabel ('\delta')


subplot(2,1,1)
rect = getrect(1);
%% getrect 从图中用鼠标截取一个矩形区域， rect 中存放的是
%% 矩形左下角的坐标 (x,y) 以及所截矩形的宽度和高度
rhomin=rect(1);　　％％截止rho
deltamin=rect(２);       ％％截止delta


％％聚类中心
%% 初始化 中心 个数
NCLUST=0;

%% cl 为归属标志数组，cl(i)=j 表示第 i 号数据点归属于第 j 个 cluster
%% 先统一将 cl 初始化为 -1
for i=1:ND
  cl(i)=-1;
end

%% 在矩形区域内统计数据点（即聚类中心）的个数
for i=1:ND
  if ( (rho(i)>rhomin) && (delta(i)>deltamin))
     NCLUST=NCLUST+1;
     cl(i)=NCLUST;                  %% 第 i 号数据点属于第 NCLUST 个 cluster
     icl(NCLUST)=i;               %% 逆映射,第 NCLUST 个 cluster 的中心为第 i 号数据点
  end
end

fprintf('NUMBER OF CLUSTERS: %i \n', NCLUST);
disp('Performing assignation')


％％将其他数据点归类
for i=1:ND
  if (cl(ordrho(i))==-1)
    cl(ordrho(i))=cl(nneigh(ordrho(i)));　　％％ordrho(i)按密度降序依次归类
  ％％nneigh(ordrho(i))密度更大距离最小的点　的类别
  end
end


%处理噪点

for i=1:ND
  halo(i)=cl(i);
end

  % 初始化数组 bord_rho 为 0,每个 cluster 定义一个 bord_rho 值
  if (NCLUST>1)
  for i=1:NCLUST
    bord_rho(i)=0.;
  end

 % 获取每一个 cluster 中平均密度的一个界 bord_rho
  for i=1:ND-1
    for j=i+1:ND
       %% 距离足够小但不属于同一个 cluster 的 i 和 j
      if ((cl(i)~=cl(j))&& (dist(i,j)<=dc))
        rho_aver=(rho(i)+rho(j))/2.;　　％％取均值密度

      ％％边界密度取两者较大值
        if (rho_aver>bord_rho(cl(i))) 
          bord_rho(cl(i))=rho_aver;　　
        end
        ％％边界密度取两者较大值
        if (rho_aver>bord_rho(cl(j))) 
          bord_rho(cl(j))=rho_aver;
        end
      end
    end
  end
  for i=1:ND
    if (rho(i)<bord_rho(cl(i)))　％％噪点
      halo(i)=0;
    end
  end
end

％％统计
for i=1:NCLUST　　
  nc=0;%% 用于累计当前 cluster 中数据点的个数
  nh=0;%% 用于累计当前 cluster 中核心数据点的个数
  for j=1:ND
    if (cl(j)==i) 
      nc=nc+1;
    end
    if (halo(j)==i) 
      nh=nh+1;
    end
  end
  fprintf('CLUSTER: %i CENTER: %i ELEMENTS: %i CORE: %i HALO: %i \n', i,icl(i),nc,nh,nc-nh);
end

cmap=colormap;
for i=1:NCLUST
   ic=int8((i*64.)/(NCLUST*1.));
   subplot(2,1,1)
   hold on
   plot(rho(icl(i)),delta(icl(i)),'o','MarkerSize',8,'MarkerFaceColor',cmap(ic,:),'MarkerEdgeColor',cmap(ic,:));
end
subplot(2,1,2)
disp('Performing 2D nonclassical multidimensional scaling')
Y1 = mdscale(dist, 2, 'criterion','metricstress');
plot(Y1(:,1),Y1(:,2),'o','MarkerSize',2,'MarkerFaceColor','k','MarkerEdgeColor','k');
title ('2D Nonclassical multidimensional scaling','FontSize',15.0)
xlabel ('X')
ylabel ('Y')
for i=1:ND
 A(i,1)=0.;
 A(i,2)=0.;
end
for i=1:NCLUST
  nn=0;
  ic=int8((i*64.)/(NCLUST*1.));
  for j=1:ND
    if (halo(j)==i)
      nn=nn+1;
      A(nn,1)=Y1(j,1);
      A(nn,2)=Y1(j,2);
    end
  end
  hold on
  plot(A(1:nn,1),A(1:nn,2),'o','MarkerSize',2,'MarkerFaceColor',cmap(ic,:),'MarkerEdgeColor',cmap(ic,:));
end

%for i=1:ND
%   if (halo(i)>0)
%      ic=int8((halo(i)*64.)/(NCLUST*1.));
%      hold on
%      plot(Y1(i,1),Y1(i,2),'o','MarkerSize',2,'MarkerFaceColor',cmap(ic,:),'MarkerEdgeColor',cmap(ic,:));
%   end
%end
faa = fopen('CLUSTER_ASSIGNATION', 'w');
disp('Generated file:CLUSTER_ASSIGNATION')
disp('column 1:element id')
disp('column 2:cluster assignation without halo control')
disp('column 3:cluster assignation with halo control')
for i=1:ND
   fprintf(faa, '%i %i %i\n',i,cl(i),halo(i));
end
