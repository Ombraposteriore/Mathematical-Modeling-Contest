clc
clear
% 文件路径
file_path = "C:\SourceCode\Python\pythonProject\BCHAIN-MKPRU.csv";  % 替换为实际的文件路径

% 读取CSV文件的B列第i到i+30行的数据

data_table = readtable(file_path);

numsi=1;
for numsi=1:7:1797
column_Value_data = data_table.Value(numsi:numsi+29);
disp("预测次数：")
disp(numsi)
%% LSTM
data=column_Value_data';  %不转置的话，无法训练lstm网络，显示维度不对。
%序列的前27个用于训练，后3个用于验证神经网络，然后往后预测7个数据。
dataTrain = data(1:30);      %定义训练集**********************************************
dataTest =  data(28:30);    %该数据是用来在最后与预测值进行对比的
%% 数据预处理
mu = mean(dataTrain);    %求均值 
sig = std(dataTrain);      %求均差 
dataTrainStandardized = (dataTrain - mu) / sig;            
%输入的每个时间步，LSTM网络学习预测下一个时间步，这里交错一个时间步效果最好。
XTrain = dataTrainStandardized(1:end-1);  
YTrain = dataTrainStandardized(2:end);  
%一维特征lstm网络训练
numFeatures = 1;   %特征为一维
numResponses = 1;  %输出也是一维
numHiddenUnits = 200;   %创建LSTM回归网络，指定LSTM层的隐含单元个数200。可调

layers = [ ...
    sequenceInputLayer(numFeatures)    %输入层
    lstmLayer(numHiddenUnits)  % lstm层，如果是构建多层的LSTM模型，可以修改。
    fullyConnectedLayer(numResponses)    %为全连接层,是输出的维数。
    regressionLayer];      %其计算回归问题的半均方误差模块 。即说明这不是在进行分类问题。

%指定训练选项，求解器设置为adam， 1000轮训练。
%梯度阈值设置为 1。指定初始学习率 0.01，在 125 轮训练后通过乘以因子 0.2 来降低学习率。
options = trainingOptions('adam', ...
    'MaxEpochs',1000, ...
    'GradientThreshold',1, ...
    'InitialLearnRate',0.01, ...      
    'LearnRateSchedule','piecewise', ...%每当经过一定数量的时期时，学习率就会乘以一个系数。
    'LearnRateDropPeriod',400, ...      %乘法之间的纪元数由“ LearnRateDropPeriod”控制。可调
    'LearnRateDropFactor',0.15, ...      %乘法因子由参“ LearnRateDropFactor”控制，可调
    'Verbose',0,  ...  %如果将其设置为true，则有关训练进度的信息将被打印到命令窗口中。默认值为true。
    'Plots','none');    %构建曲线图 将'training-progress'替换为none
net = trainNetwork(XTrain,YTrain,layers,options); 
%% 神经网络初始化
net = predictAndUpdateState(net,XTrain);  %将新的XTrain数据用在网络上进行初始化网络状态
[net,YPred] = predictAndUpdateState(net,YTrain(end));  %用训练的最后一步来进行预测第一个预测值，给定一个初始值。这是用预测值更新网络状态特有的。
%% 进行用于验证神经网络的数据预测（用预测值更新网络状态）
for i = 2:7  %从第二步开始，这里进行7次单步预测
    [net,YPred(:,i)] = predictAndUpdateState(net,YPred(:,i-1),'ExecutionEnvironment','cpu');  %predictAndUpdateState函数是一次预测一个值并更新网络状态
end
%% 验证神经网络
YPred = sig*YPred + mu;      %使用先前计算的参数对预测去标准化。
rmse = sqrt(mean((YPred(1:3)-dataTest).^2)) ;     %计算均方根误差 (RMSE)。
result=YPred(1:7)';

%% 写入forecast列
data_table.forecast(numsi+30:numsi+36) = result;

% 将更新后的数据写回文件
writetable(data_table, file_path);

%% free
clear YPred
clear YTrian
clear XTrain
clear YTrain
clear sig
clear rmse
clear result
clear options
clear numResponses

clear numFeatures
clear numHiddenUnits
clear mu
clear net
clear layers
clear i
clear dataTrainStandardized
clear dataTrain
clear dataTest
clear data
clear column_Value_data
 end
