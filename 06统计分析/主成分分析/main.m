clear all;
clc
load data1.mat   % 主成分聚类
% load data2.mat   % 主成分回归

% 注意，这里可以对数据先进行描述性统计
% 描述性统计的内容见第5讲.相关系数
[n,p] = size(x);   % n是样本个数，p是指标个数

%% 第1步：对数据x进行标准化处理
%std (x, flag,dim)
% fla表示标注公差时是要除以n还是n－1
% flag==0.........是除以n－1
% flag==1.........是除以n
% dim表示维数
% dim==1..........是按照列分
% dim==2..........是按照行分 若是三维的矩阵，dim=＝3就按照第三维来分数据
% 默认std格式是std(x,0,1);

% 标准化公式为 (每列数据的均值-样本标准差)/标准差
X=zscore(x);       % matlab内置的标准化函数（x-mean(x)）/std(x)

%% 第2步：计算样本协方差矩阵
% R = cov(X);      % 求协方差矩阵

%% 注意：以上两步可合并为下面一步：直接计算样本相关系数矩阵
R = corrcoef(x);
disp('样本相关系数矩阵为：') 
disp(R)

%% 第三步：计算R的特征值和特征向量
% 注意：R是半正定矩阵，所以其特征值不为负数
% R同时是对称矩阵，Matlab计算对称矩阵时，会将特征值按照从小到大排列哦
% eig函数的详解见第一讲层次分析法的视频
[V,D] = eig(R);   % V 特征向量矩阵  D 特征值构成的对角矩阵


%% 第四步：计算主成分贡献率和累计贡献率
PV = diag(D);       % diag函数用于得到一个矩阵的主对角线元素值(返回的是列向量)  PV ：flag value 特征值
PV = PV(end:-1:1);  % 因为特征向量是从小大到排序的，为了方便计算，翻转一下
contribution_rate = PV / sum(PV);              % 计算贡献率
cum_contribution_rate = cumsum(PV)/ sum(PV);   % 计算累计贡献率  cumsum是求累加值的函数

disp('特征值为：')
disp(PV')           %转置为行向量，方便展示
disp('贡献率为：')
disp(contribution_rate')
disp('累计贡献率为：')
disp(cum_contribution_rate')
disp('与特征值对应的特征向量矩阵为：')
%  注意：这里的特征向量要和特征值一一对应，之前特征值相当于颠倒过来了，因此特征向量的各列需要颠倒过来
%  rot90函数可以使一个矩阵逆时针旋转90度，然后再转置，就可以实现将矩阵的列颠倒的效果
V=rot90(V)';
disp(V)


%% 计算我们所需要的主成分的值
m = input('请输入需要保存的主成分的个数:  ');
F = zeros(n,m);  %初始化保存主成分的矩阵（每一列是一个主成分）
for i = 1:m
    ai = V(:,i)';               % 将第i个特征向量取出，并转置为行向量
    Ai = repmat(ai,n,1);        % 将这个行向量重复n次，构成一个n*p的矩阵
    F(:, i) = sum(Ai .* X, 2);  % 注意，对标准化的数据求了权重后要计算每一行的和
end

% (1)主成分聚类 ： 将主成分指标所在的F矩阵复制到Excel表格，然后再用Spss进行聚类
% 在Excel第一行输入指标名称（F1,F2, ..., Fm）
% 双击Matlab工作区的F,进入变量编辑中，然后复制里面的数据到Excel表格
% 导出数据之后，我们后续的分析就可以在Spss中进行。
% 
% （2）主成分回归：将x使用主成分得到主成分指标，并将y标准化，接着导出到Excel，然后再使用Stata回归
% Y = zscore(y);  % 一定要将y进行标准化哦~
% 在Excel第一行输入指标名称（Y,F1, F2, ..., Fm）
% 分别双击Matlab工作区的Y和F,进入变量编辑中，然后复制里面的数据到Excel表格
% 导出数据之后，我们后续的分析就可以在Stata中进行。






