clear all;
clc
load  data1
plot(x,y,'o')    %画样本
% 给x和y轴加上标签
xlabel('x的值')
ylabel('y的值')
n = size(x,1);
%  y=kxd+b

k = (n*sum(x.*y)-sum(x)*sum(y))/(n*sum(x.*x)-sum(x)*sum(x)) 
b = (sum(x.*x)*sum(y)-sum(x)*sum(x.*y))/(n*sum(x.*x)-sum(x)*sum(x))
hold on %继续在之前的图形上来画图形
grid on %显示网格线

% % 画出y=kx+b的函数图像 plot(x,y)
% % 传统的画法：模拟生成x和y的序列，比如要画出[0,5]上的图形
% xx = 2.5: 0.1 :7  % 间隔设置的越小画出来的图形越准确
% yy = k * xx + b  % k和b都是已知值
% plot(xx,yy,'-')

% 匿名函数的基本用法。
% handle = @(arglist) anonymous_function
% 其中handle为调用匿名函数时使用的名字。
% arglist为匿名函数的输入参数，可以是一个，也可以是多个，用逗号分隔。
% anonymous_function为匿名函数的表达式。
% 举个小例子
%  z=@(x,y) x^2+y^2; 
%  z(1,2) 
% % ans =  5
% fplot函数可用于画出匿名一元函数的图形。
% fplot(f,xinterval) 将匿名函数f在指定区间xinterval绘图。xinterval =  [xmin xmax] 表示定义域的范围

f=@(x) k*x+b;
fplot(f,[2.5,7]);
legend('样本数据','拟合函数','location','SouthEast')

y_hat = k*x+b; %y的拟合值
% 只有是线性函数（线性于参数）才能用R^2（拟合优度，R^2越接近1，说明拟合度越好）
% 复杂的函数比如指数函数，直接看SSE即可

SSR = sum((y_hat-mean(y)).^2)  %回归平方和 mean()是求均值的函数
SSE = sum((y_hat-y).^2)        %误差平方和
SST = sum((y-mean(y)).^2)      %总体平方和
SST-SSE-SSR   %5.6843e-14=5.6843*10^-14  matlab浮点数计算的一个误差，就是等于零
R_2 = SSR / SST



