clear;
clc;
close all
tStart = tic; % 算法计时器

%% 自定义参数

data= load('citys_data.mat');
cityNum = size(data.citys,1);
cities = data.citys';
%cityNum = 100;
maxGEN = 1000;
popSize = 100; % 遗传算法种群大小
crossoverProbabilty = 0.9; %交叉概率
mutationProbabilty = 0.1; %变异概率
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
gbest = Inf;
% 随机生成城市位置
%cities = rand(2,cityNum) * 100;%100是最远距离
% 计算上述生成的城市距离
distances = calculateDistance(cities);
% 生成种群，每个个体代表一个路径
pop = zeros(popSize, cityNum);
for i=1:popSize
    pop(i,:) = randperm(cityNum);
end
offspring = zeros(popSize,cityNum);
%保存每代的最小路径便于画图
minPathes = zeros(maxGEN,1);
% GA算法
for  gen=1:maxGEN
    % 计算适应度的值，即路径总距离
    [fval, sumDistance, minPath, maxPath] = fitness(distances, pop);
    % 轮盘赌选择
    tournamentSize=4; %设置大小
    for k=1:popSize
        % 选择父代进行交叉
        tourPopDistances=zeros( tournamentSize,1);
        for i=1:tournamentSize
            randomRow = randi(popSize);
            tourPopDistances(i,1) = sumDistance(randomRow,1);
        end
        % 选择最好的，即距离最小的
        parent1  = min(tourPopDistances);
        [parent1X,parent1Y] = find(sumDistance==parent1,1, 'first');
        parent1Path = pop(parent1X(1,1),:);
        for i=1:tournamentSize
            randomRow = randi(popSize);
            tourPopDistances(i,1) = sumDistance(randomRow,1);
        end
        parent2  = min(tourPopDistances);
        [parent2X,parent2Y] = find(sumDistance==parent2,1, 'first');
        parent2Path = pop(parent2X(1,1),:);
        subPath = crossover(parent1Path, parent2Path, crossoverProbabilty);%交叉
        subPath = mutate(subPath, mutationProbabilty);%变异
        offspring(k,:) = subPath(1,:);
        minPathes(gen,1) = minPath;
    end
    fprintf('代数:%d   最短路径:%.2fKM \n', gen,minPath);
    % 更新
    pop = offspring;
    % 画出当前状态下的最短路径
    if minPath < gbest
        gbest = minPath;
        paint(cities, pop, gbest, sumDistance,gen);
    end
end
figure 
plot(minPathes, 'MarkerFaceColor', 'red','LineWidth',1);
title('收敛曲线图（每一代的最短路径）');
set(gca,'ytick',500:100:5000); 
ylabel('路径长度');
xlabel('迭代次数');
grid on
tEnd = toc(tStart);
fprintf('时间:%d 分  %f 秒.\n', floor(tEnd/60), rem(tEnd,60));

function [ distances ] = calculateDistance( city )
%计算距离
 
    [~, col] = size(city);
    distances = zeros(col);
    for i=1:col
        for j=1:col
            distances(i,j)= distances(i,j)+ sqrt( (city(1,i)-city(1,j))^2 + (city(2,i)-city(2,j))^2  );           
        end
    end
end

function [childPath] = crossover(parent1Path, parent2Path, prob)
% 交叉
 
    random = rand();
    if prob >= random
        [l, length] = size(parent1Path);
        childPath = zeros(l,length);
        setSize = floor(length/2) -1;
        offset = randi(setSize);
        for i=offset:setSize+offset-1
            childPath(1,i) = parent1Path(1,i);
        end
        iterator = i+1;
        j = iterator;
        while any(childPath == 0)
            if j > length
                j = 1;
            end
 
            if iterator > length
                iterator = 1;
            end
            if ~any(childPath == parent2Path(1,j))
                childPath(1,iterator) = parent2Path(1,j);
               iterator = iterator + 1;
            end
            j = j + 1;
        end
    else
        childPath = parent1Path;
    end
end

function [ fitnessvar, sumDistances,minPath, maxPath ] = fitness( distances, pop )
% 计算整个种群的适应度值
 
    [popSize, col] = size(pop);
    sumDistances = zeros(popSize,1);
    fitnessvar = zeros(popSize,1);
    for i=1:popSize
       for j=1:col-1
          sumDistances(i) = sumDistances(i) + distances(pop(i,j),pop(i,j+1));
       end 
    end
    minPath = min(sumDistances);
    maxPath = max(sumDistances);
    for i=1:length(sumDistances)
        fitnessvar(i,1)=(maxPath - sumDistances(i,1)+0.000001) / (maxPath-minPath+0.00000001);
    end
end

%交叉
function [ mutatedPath ] = mutate( path, prob )
%对指定的路径利用指定的概率进行更新
 
    random = rand();
    if random <= prob
        [l,length] = size(path);
        index1 = randi(length);
        index2 = randi(length);
        %交换
        temp = path(l,index1);
        path(l,index1) = path(l,index2);
        path(l,index2)=temp;
    end
        mutatedPath = path; 
end

%paint
function [ output_args ] = paint( cities, pop, minPath, totalDistances,gen)
    gNumber=gen;
    [~, length] = size(cities);
    xDots = cities(1,:);
    yDots = cities(2,:);
    %figure(1);
    title('GA TSP');
    plot(xDots,yDots, 'p', 'MarkerSize', 14, 'MarkerFaceColor', 'blue');
    xlabel('经度');
    ylabel('纬度');
    axis equal
    hold on
    [minPathX,~] = find(totalDistances==minPath,1, 'first');
    bestPopPath = pop(minPathX, :);
    bestX = zeros(1,length);
    bestY = zeros(1,length);
    for j=1:length
       bestX(1,j) = cities(1,bestPopPath(1,j));
       bestY(1,j) = cities(2,bestPopPath(1,j));
    end
    title('GA TSP');
    plot(bestX(1,:),bestY(1,:), 'red', 'LineWidth', 1.25);
    legend('城市', '路径');
    axis equal
    grid on
    %text(5,0,sprintf('迭代次数: %d 总路径长度: %.2f',gNumber, minPath),'FontSize',10);
    drawnow
    hold off
end