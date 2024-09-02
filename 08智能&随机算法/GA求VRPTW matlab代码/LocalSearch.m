
%% 局部搜索函数
%输入：SelCh               被选择的个体
%输入：cusnum              顾客数目
%输入：cap                 最大载重量
%输入：demands             需求量
%输入：a                   顾客时间窗开始时间[a[i],b[i]]
%输入：b                   顾客时间窗结束时间[a[i],b[i]]
%输入：L                   配送中心时间窗结束时间
%输入：s                   客户点的服务时间
%输入：dist                距离矩阵，满足三角关系，暂用距离表示花费c[i][j]=dist[i][j]
%输出：SelCh               进化逆转后的个体
function SelCh=LocalSearch(SelCh,cusnum,cap,demands,a,b,L,s,dist,alpha,belta)
D=15;                                                       %Remove过程中的随机元素
toRemove=15;                                                %将要移出顾客的数量
[row,N]=size(SelCh);
for i=1:row
    [VC,NV,TD,violate_num,violate_cus]=decode(SelCh(i,:),cusnum,cap,demands,a,b,L,s,dist);
    CF=costFuction(VC,a,b,s,L,dist,demands,cap,alpha,belta);
    %Remove
    [removed,rfvc] = Remove(cusnum,toRemove,D,dist,VC);
    %Re-inserting
    [ReIfvc,RTD] = Re_inserting(removed,rfvc,L,a,b,s,dist,demands,cap);
    %计算惩罚函数
    RCF=costFuction(ReIfvc,a,b,s,L,dist,demands,cap,alpha,belta);
    if RCF<CF
        chrom=change(ReIfvc,N,cusnum);
        if length(chrom)~=N
            record=ReIfvc;
            break
        end
        SelCh(i,:)=chrom;
    end
end
