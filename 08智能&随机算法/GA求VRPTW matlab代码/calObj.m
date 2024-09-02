
%% 计算种群的目标函数值
%输入：Chrom               种群
%输入：cusnum              顾客数目
%输入：cap                 最大载重量
%输入：demands             需求量
%输入：a                   顾客时间窗开始时间[a[i],b[i]]
%输入：b                   顾客时间窗结束时间[a[i],b[i]]
%输入：L                   配送中心时间窗结束时间
%输入：s                   客户点的服务时间
%输入：dist                距离矩阵，满足三角关系，暂用距离表示花费c[i][j]=dist[i][j]
%输出：ObjV                每个个体的目标函数值，定义为车辆使用数目*10000+车辆行驶总距离
function ObjV=calObj(Chrom,cusnum,cap,demands,a,b,L,s,dist,alpha,belta)
% route=1:cusnum;
% G= part_length(route,dist);               %G为对每条不可行路径的惩罚权重
NIND=size(Chrom,1);                         %种群数目
ObjV=zeros(NIND,1);                         %储存每个个体函数值
G=10;                                       %G为对每条不可行路径的惩罚权重
for i=1:NIND
    [VC,NV,TD,violate_num,violate_cus]=decode(Chrom(i,:),cusnum,cap,demands,a,b,L,s,dist);
    costF=costFuction(VC,a,b,s,L,dist,demands,cap,alpha,belta);
    ObjV(i)=costF;
%     ObjV(i)=NV+G*violate_cus;
end
end

