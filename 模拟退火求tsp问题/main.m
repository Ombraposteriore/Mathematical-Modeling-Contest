%% 模拟退火解决TSP问题
%main.m主函数
tic
clear;clc
coord = [26.99611373,17.4952427;30.86325116,42.42106760;43.61357686,38.9299278;
    35.00118717,36.98696184;20.20818722,46.55925669;46.18764943,42.83843589;
    40.6920241,31.94898546;32.50662006,10.62166377;45.88614135,13.7874543;
    24.72317638,30.59033166;23.50924028,48.09930553;34.75421424,43.72244785;
    13.4758176,47.40048984;25.01990127,48.9365234];%坐标信息自己改

n = size(coord,1);  % 城市的数目
 
figure  % 新建一个图形窗口
plot(coord(:,1),coord(:,2),'o');   % 画出城市的分布散点图
hold on % 等一下要接着在这个图形上画图的
 
d = zeros(n);   % 初始化两个城市的距离矩阵全为0
for i = 2:n  
    for j = 1:i  
        coord_i = coord(i,:);   x_i = coord_i(1);     y_i = coord_i(2);  % 城市i的横坐标为x_i，纵坐标为y_i
        coord_j = coord(j,:);   x_j = coord_j(1);     y_j = coord_j(2);  % 城市j的横坐标为x_j，纵坐标为y_j
        d(i,j) = sqrt((x_i-x_j)^2 + (y_i-y_j)^2);   % 计算城市i和j的距离
    end
end
d = d+d';   % 生成距离矩阵的对称的一面
 
%% 参数初始化
T0 = 1000;   % 初始温度
T = T0; % 迭代中温度会发生改变，第一次迭代时温度就是T0
maxgen = 1000;  % 最大迭代次数
Lk = 500;  % 每个温度下的迭代次数
alpfa = 0.95;  % 温度衰减系数
 
%%  随机生成一个初始解
path0 = randperm(n);  % 生成一个1-n的随机打乱的序列作为初始的路径
result0 = calculate_tsp_d(path0,d); % 调用我们自己写的calculate_tsp_d函数计算当前路径的距离
 
%% 定义一些保存中间过程的量，方便输出结果和画图
min_result = result0;     % 初始化找到的最佳的解对应的距离为result0
RESULT = zeros(maxgen,1); % 记录每一次外层循环结束后找到的min_result (方便画图）
 
%% 模拟退火过程
for iter = 1 : maxgen  % 外循环, 我这里采用的是指定最大迭代次数
    for i = 1 : Lk  %  内循环，在每个温度下开始迭代
        path1 = gen_new_path(path0);  % 调用我们自己写的gen_new_path函数生成新的路径
        result1 = calculate_tsp_d(path1,d); % 计算新路径的距离
        %如果新解距离短，则直接把新解的值赋值给原来的解
        if result1 < result0    
            path0 = path1; % 更新当前路径为新路径
            result0 = result1; 
        else
            p = exp(-(result1 - result0)/T); % 根据Metropolis准则计算一个概率
            if rand(1) < p   % 生成一个随机数和这个概率比较，如果该随机数小于这个概率
                path0 = path1;  % 更新当前路径为新路径
                result0 = result1; 
            end
        end
        % 判断是否要更新找到的最佳的解
        if result0 < min_result  % 如果当前解更好，则对其进行更新
            min_result = result0;  % 更新最小的距离
            best_path = path0;  % 更新找到的最短路径
        end
    end
    RESULT(iter) = min_result; % 保存本轮外循环结束后找到的最小距离
    T = alpfa*T;   % 温度下降       
end
 
disp('最佳的方案是：'); disp(mat2str(best_path))
disp('此时最优值是：'); disp(min_result)
 
best_path = [best_path,best_path(1)];   % 在最短路径的最后面加上一个元素，即第一个点（我们要生成一个封闭的图形）
n = n+1;  % 城市的个数加一个（紧随着上一步）
for i = 1:n-1 
    j = i+1;
    coord_i = coord(best_path(i),:);   x_i = coord_i(1);     y_i = coord_i(2); 
    coord_j = coord(best_path(j),:);   x_j = coord_j(1);     y_j = coord_j(2);
    plot([x_i,x_j],[y_i,y_j],'-b')    % 每两个点就作出一条线段，直到所有的城市都走完
%     pause(0.02)  % 暂停0.02s再画下一条线段
    hold on
end
 
%% 画出每次迭代后找到的最短路径的图形
figure
plot(1:maxgen,RESULT,'b-');
xlabel('迭代次数');
ylabel('最短路径');
toc
