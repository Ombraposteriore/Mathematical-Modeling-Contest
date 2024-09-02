
clear
clc
close all
tic
%% 用importdata这个函数来读取文件
c101=importdata('c101.txt');
cap=500;                                                        %车辆最大装载量
%% 提取数据信息
E=c101(1,5);                                                    %配送中心时间窗开始时间
L=c101(1,6);                                                    %配送中心时间窗结束时间
vertexs=c101(:,2:3);                                            %所有点的坐标x和y
customer=vertexs(2:end,:);                                       %顾客坐标
cusnum=size(customer,1);                                         %顾客数
v_num=2;                                                        %车辆最多使用数目
demands=c101(2:end,4);                                          %需求量
a=c101(2:end,5);                                                %顾客时间窗开始时间[a[i],b[i]]
b=c101(2:end,6);                                                %顾客时间窗结束时间[a[i],b[i]]
s=c101(2:end,7);                                                %客户点的服务时间
% 添加所有的节点
nodes = 16;
% 添加所有的边
edges=[1    1	0	0	1	0	1	0	0	0	1	0	0	0	0	0;
1	1	1	0	1	0	0	0	0	0	0	0	0	0	0	0;
0	1	1	1	1	1	0	0	0	0	0	0	0	0	0	0;
0	0	1	1	0	1	0	0	0	1	0	0	0	0	0	0;
1	1	1	0	1	1	1	1	1	0	0	0	0	0	0	0;
0	0	1	1	1	1	0	0	1	1	0	0	0	0	0	0;
1	0	0	0	1	0	1	1	0	0	1	1	0	0	0	0;
0	0	0	0	1	0	1	1	1	0	0	1	1	0	1	0;
0	0	0	0	1	1	0	1	1	1	0	1	1	1	1	1;
0	0	0	1	0	1	0	0	1	1	0	0	0	1	0	1;
1	0	0	0	0	0	1	0	0	0	1	1	0	0	0	0;
0	0	0	0	0	0	1	1	1	0	1	1	1	0	1	0;
0	0	0	0	0	0	0	1	1	0	0	1	1	1	1	1;
0	0	0	0	0	0	0	0	1	1	0	0	1	1	1	1;
0	0	0	0	0	0	0	1	1	0	0	1	1	1	1	0;
0	0	0	0	0	0	0	0	1	1	0	0	1	1	0	1
];
% 客户位置和需求量矩阵 ，默认第一个为仓库
locations =[20.20818722,46.55925669;30.86325116,42.42106760;43.61357686,38.9299278;
    35.00118717,36.98696184;26.99611373,17.4952427;46.18764943,42.83843589;
    40.6920241,31.94898546;32.50662006,10.62166377;45.88614135,13.7874543;
    24.72317638,30.59033166;23.50924028,48.09930553;34.75421424,43.72244785;
    13.4758176,47.40048984;25.01990127,48.9365234;13.4758176,47.40048984;25.01990127,48.9365234];
%   计算距离矩阵
dist = zeros(nodes, nodes);
for i = 1:nodes
    for j = 1:nodes
        if edges(i,j)==1  % 判断是否有直接连接
            dist(i,j) = norm(locations(i,:) - locations(j,:));
        else
            dist(i,j) = 100000;  % 如果不能直接连接，将距离设为无限大
        end
    end
end 
% h=pdist(vertexs);
% dist=squareform(h);                                             %距离矩阵，满足三角关系，暂用距离表示花费c[i][j]=dist[i][j]
%% 遗传算法参数设置
alpha=1000;                                                       %违反的容量约束的惩罚函数系数
belta=0;                                                      %违反时间窗约束的惩罚函数系数
NIND=300;                                                       %种群大小
MAXGEN=800;                                                     %迭代次数
Pc=0.91;                                                         %交叉概率
Pm=0.01;                                                        %变异概率
GGAP=0.9;                                                       %代沟(Generation gap)
N=cusnum+v_num-1;                                               %染色体长度=顾客数目+车辆最多使用数目-1
%% 初始化种群
init_vc=init(cusnum,a,demands,cap);                             %构造初始解
Chrom=InitPopCW(NIND,N,cusnum,init_vc);
%% 输出随机解的路线和总距离
disp('初始种群中的一个随机值:')
[VC,NV,TD,violate_num,violate_cus]=decode(Chrom(1,:),cusnum,cap,demands,a,b,L,s,dist);
% disp(['总距离：',num2str(TD)]);
disp(['车辆使用数目：',num2str(NV),'，车辆行驶总距离：',num2str(TD),'，违反约束路径数目：',num2str(violate_num),'，违反约束顾客数目：',num2str(violate_cus)]);
disp('~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~')
%% 优化
gen=1;
figure;
hold on;box on
xlim([0,MAXGEN])
title('优化过程')
xlabel('代数')
ylabel('最优值')
ObjV=calObj(Chrom,cusnum,cap,demands,a,b,L,s,dist,alpha,belta);             %计算种群目标函数值
preObjV=min(ObjV);
while gen<=MAXGEN
    %% 计算适应度
    ObjV=calObj(Chrom,cusnum,cap,demands,a,b,L,s,dist,alpha,belta);             %计算种群目标函数值
    line([gen-1,gen],[preObjV,min(ObjV)]);pause(0.0001)
    preObjV=min(ObjV);
    FitnV=Fitness(ObjV);
    %% 选择
    SelCh=Select(Chrom,FitnV,GGAP);
    %% OX交叉操作
    SelCh=Recombin(SelCh,Pc);
    %% 变异
    SelCh=Mutate(SelCh,Pm);
    %% 局部搜索操作
    SelCh=LocalSearch(SelCh,cusnum,cap,demands,a,b,L,s,dist,alpha,belta);
    %% 重插入子代的新种群
    Chrom=Reins(Chrom,SelCh,ObjV);
    %% 删除种群中重复个体，并补齐删除的个体
    Chrom=deal_Repeat(Chrom);
    %% 打印当前最优解
    ObjV=calObj(Chrom,cusnum,cap,demands,a,b,L,s,dist,alpha,belta);             %计算种群目标函数值
    [minObjV,minInd]=min(ObjV);
    disp(['第',num2str(gen),'代最优解:'])
    [bestVC,bestNV,bestTD,best_vionum,best_viocus]=decode(Chrom(minInd(1),:),cusnum,cap,demands,a,b,L,s,dist);
    disp(['车辆使用数目：',num2str(bestNV),'，车辆行驶总距离：',num2str(bestTD),'，违反约束路径数目：',num2str(best_vionum),'，违反约束顾客数目：',num2str(best_viocus)]);
    fprintf('\n')
    %% 更新迭代次数
    gen=gen+1 ;
end
%% 画出最优解的路线图
ObjV=calObj(Chrom,cusnum,cap,demands,a,b,L,s,dist,alpha,belta);             %计算种群目标函数值
[minObjV,minInd]=min(ObjV);
%% 输出最优解的路线和总距离
disp('最优解:')
bestChrom=Chrom(minInd(1),:);
[bestVC,bestNV,bestTD,best_vionum,best_viocus]=decode(bestChrom,cusnum,cap,demands,a,b,L,s,dist);
disp(['车辆使用数目：',num2str(bestNV),'，车辆行驶总距离：',num2str(bestTD),'，违反约束路径数目：',num2str(best_vionum),'，违反约束顾客数目：',num2str(best_viocus)]);
disp('-------------------------------------------------------------')
%% 判断最优解是否满足时间窗约束和载重量约束，0表示违反约束，1表示满足全部约束
flag=Judge(bestVC,cap,demands,a,b,L,s,dist);
%% 检查最优解中是否存在元素丢失的情况，丢失元素，如果没有则为空
DEL=Judge_Del(bestVC);
%% 画出最终路线图
draw_Best(bestVC,vertexs);
save c101.mat
toc