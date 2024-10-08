% 1.输入判断矩阵
% 2.特征值法求权重(将最大的特征向量作为权向量)
% 3.一致性检验

clear all
clc
disp('请输入判断矩阵A');
A=input('A=');
[n,n]=size(A);%先求出A的阶数，后面算一致性检验需要


[V,D] = eig(A);       %V是特征向量，D是由特征值构成对角矩阵，V的第1列特征向量对应D的第1列的特征值
disp('特征向量');
disp(V);
disp('特征值');
disp(max(D));

DMax= max(max(D));     %求出最大特征值 也可以写成DMax= max(D(:));
%既然找到了最大特征值，怎么找对应特征向量？
%1.比较 D == DMax，返回一个只含bool值的矩阵(0,1)  D上的元素只有等于DMax(最大特征值)上的位置为1，其余位置为零
%2.利用find函数找到最大特征值在D中的索引
[r,c]=find(D == DMax , 1);
disp('最大的特征值DMax在D中的列数为：');
disp(c);
disp('最大特征值DMax对应的特征向量为：');
disp(V(:,c));

disp('最终权向量为：');
disp( V(:,c) ./ sum(V(:,c)) )  %特征向量的归一化


%计算一致性比例CR的环节
CI = (DMax - n) / (n-1);
RI=[0 0.0001 0.52 0.89 1.12 1.26 1.36 1.41 1.46 1.49 1.52 1.54 1.56 1.58 1.59];  %这里的RI最多支持 n = 15
%n=2时，2阶一定是一致性矩阵，为了避免分母为零，取第二个数为比零大一点的数
CR=CI/RI(n);
disp('一致性指标CI=');disp(CI);
disp('一致性比例CR=');disp(CR);
if CR<0.10
    disp('因为CR<0.10，A矩阵通过一致性检验。');
else
    disp('CR >= 0.10，判断矩阵A需要进行修改!');
end


% A=[1    1/2   4   3    3
%    2    1     7   5    5
%    1/4  1/7   1   1/2  1/3
%    1/3  1/5   2   1     1
%    1/3  1/5   3   1     1]
