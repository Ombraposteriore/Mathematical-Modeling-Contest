%不采用GA工具，手撸代码求解一元函数的最大值
%hzj
%2023-08-11


%标准遗传算法求函数极值
%% 初始化参数
clear all; %清除所有变量
close all; %清图
clc; %清屏
NP = 50; %种群数量
L = 20; %二进制位串长度
Pc = 0.8; %交叉率
Pm = 0.1; %变异率
G = 100; %最大遗传代数
Xs = 10; %上限
Xx = 0; %下限
f = rand(NP,L); %随机获得初始种群


%%%%%%%%%%%%%%%%%%遗传算法循环%%%%%%%%%%%%%%%%
for k = 1:G
    %%%%%%%%%%将二进制解码为定义域范围内十进制%%%%%%%%%%
    for i = 1:NP
        U = f(i,:);%第i行所有元素
        m = 0;
        for j = 1:L
            m = U(j)*2^(j-1)+m;
        end
        x(i) = Xx+m*(Xs-Xx)/(2^L-1);
        Fit(i) = func1(x(i));
    end
    maxFit = max(Fit); %最大值
    minFit = min(Fit); %最小值
    rr = find(Fit==maxFit);
    fBest = f(rr(1,1),:); %历代最优个体
    xBest = x(rr(1,1));
    Fit = (Fit-minFit)/(maxFit-minFit); %归一化适应度值

    %======轮赌盘========
    %%%%%%%%%%%%%%基于轮盘赌的复制操作%%%%%%%%%%%%%
    sum_Fit = sum(Fit);
    fitvalue = Fit./sum_Fit;
    fitvalue = cumsum(fitvalue);
    ms = sort(rand(NP,1));
    fiti = 1;
    newi = 1;
    while newi <= NP
        if (ms(newi)) < fitvalue(fiti)
            nf(newi,:) = f(fiti,:);
            newi = newi+1;
        else
            fiti = fiti+1;
        end
    end

      %%%%%%%%%%%%%%%基于概率的交叉操作%%%%%%%%%%%%%
    for i = 1:2:NP
        p = rand;
        if p < Pc
            q = rand(1,L);
            for j = 1:L
                if q(j)==1;
                    temp = nf(i+1,j);
                    nf(i+1,j) = nf(i,j);
                    nf(i,j) = temp;
                end
            end
        end
    end

     %%%%%%%%%%%%%基于概率的变异操作%%%%%%%%%%%%%%
    i = 1;
    while i <= round(NP*Pc)
        h = randi([1,NP],1,1); %随机选取一个需要变异的染色体
        for j = 1:round(L*Pc)
            g = randi([1,L],1,1); %随机选取需要变异的基因数
            nf(h,g) =~ nf(h,g);
        end
        i = i+1;
    end
    f = nf;
    f(1,:) = fBest; %保留最优个体在新种群中
    trace(k) = maxFit; %历代最优适应度
end
xBest; %最优个体
figure
plot(trace)
xlabel('迭代次数')
ylabel('目标函数值')
title('适应度进化曲线')

x=0:0.1:10;
y=x+10*sin(5*x)+7*cos(4*x);
plot(x,y);hold on
yBest = xBest+10*sin(5*xBest)+7*cos(4*xBest);
plot(xBest,yBest,'r*');
hold off

%% 适应度函数
function result = func1(x)
fit = x+10*sin(5*x)+7*cos(4*x);
result = fit;
end