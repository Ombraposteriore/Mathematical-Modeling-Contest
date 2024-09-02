%1991年各省消费数据，采用libsvm分类
%2023-8-9
%hzj

clc;
clear;
close all

%数据导入
load("data_1991.mat")
%数据标准化
matrix_std = mapstd(matrix);

%数据划分
matrix_train = matrix_std(1:27,:);
label_train = label(1:27);
matrix_test = matrix_std(28:end,:);
label_test = label(28:end);

%训练
[c,g] = meshgrid(-10:0.2:10,-10:0.2:10);
[m,n] = size(c);
cg = zeros(m,n);
eps = 10^(-4);
v = 5;
bestc = 1;
bestg = 0.1;
bestacc = 0;
for i = 1:m
    for j = 1:n
        cmd = ['-v ',num2str(v),' -t 2',' -c ',num2str(2^c(i,j)),' -g ',num2str(2^g(i,j))];
        cg(i,j) = libsvmtrain(label_train,matrix_train,cmd);    
        if cg(i,j) > bestacc
            bestacc = cg(i,j);
            bestc = 2^c(i,j);
            bestg = 2^g(i,j);
        end       
        if abs( cg(i,j)-bestacc )<=eps && bestc > 2^c(i,j)
            bestacc = cg(i,j);
            bestc = 2^c(i,j);
            bestg = 2^g(i,j);
        end              
    end
end
cmd = [' -t 2',' -c ',num2str(bestc),' -g ',num2str(bestg)];

% model = libsvmtrain(label_train,matrix_train);
% model = libsvmtrain(label_train,matrix_train,'-s 1 -t 2 -c 0.2');%对比调参后的结果
model = libsvmtrain(label_train,matrix_train,cmd);%对比调参后的结果 
%预测
%[predict_label,acc,pro]=libsvmpredict(label_test,matrix_test,model);

%预测
matrix_t= matrix_std(15:26,:);
label_t = label(15:26);
[predict_label,acc,pro]=libsvmpredict(label_t,matrix_t,model);




