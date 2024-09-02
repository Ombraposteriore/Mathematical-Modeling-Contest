
%% 初始化种群
%输入NIND：                种群大小
%输入N：                   染色体长度
%输入：cusnum              顾客数目
%输入init_vc：             初始配送方案
%输出Chrom：               初始种群

function Chrom=InitPopCW(NIND,N,cusnum,init_vc)
Chrom=zeros(NIND,N);%用于存储种群
chrom=change(init_vc,N,cusnum);
for j=1:NIND
    Chrom(j,:)=chrom;
end