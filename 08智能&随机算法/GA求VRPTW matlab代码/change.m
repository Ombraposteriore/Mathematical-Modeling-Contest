
%% 配送方案与个体之间进行转换
function chrom=change(VC,N,cusnum)
NV=size(VC,1);                         %车辆使用数目
chrom=[];
for i=1:NV
    if (cusnum+i)<=N
        chrom=[chrom,VC{i},cusnum+i];
    else
        chrom=[chrom,VC{i}];
    end
end
if length(chrom)<N                          %如果染色体长度小于N，则需要向染色体添加配送中心编号
    supply=(cusnum+NV+1):N;
    chrom=[chrom,supply];
end

end

