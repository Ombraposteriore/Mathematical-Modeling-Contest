
% ����SVM�Ļع�Ԥ�����������ָ֤������ָ��Ԥ��
% by hzj
% 2023.08.10
%% ��ջ�������
function chapter_sh
tic;         %��ʼ��ʱ
close all;
%clear;
clc;
format compact;  %ѹ����ʽ�����ض���Ļ��з�
%% ���ݵ���ȡ��Ԥ����

% �������������ָ֤��(1990.12.19-2009.08.19)
% ������һ��4579*6��double�͵ľ���,ÿһ�б�ʾÿһ�����ָ֤��
% 6�зֱ��ʾ������ָ֤���Ŀ���ָ��,ָ�����ֵ,ָ�����ֵ,����ָ��,���ս�����,���ս��׶�.
load chapter_sh.mat;

% ��ȡ����
[m,n] = size(sh);%m=4579,n=6
ts = sh(2:m,1);
tsx = sh(1:m-1,:);

% ����ԭʼ��ָ֤����ÿ�տ�����
figure;
plot(ts,'LineWidth',2);
title('��ָ֤����ÿ�տ�����(1990.12.20-2009.08.19)','FontSize',12);
xlabel('����������(1990.12.19-2009.08.19)','FontSize',12);
ylabel('������','FontSize',12);
grid on;

% ����Ԥ����,��ԭʼ���ݽ��й�һ��
ts = ts';
tsx = tsx';

% mapminmaxΪmatlab�Դ���ӳ�亯��	
% ��ts���й�һ��
[TS,TSps] = mapminmax(ts,1,2);	

% ����ԭʼ��ָ֤����ÿ�տ�������һ�����ͼ��
figure;
plot(TS,'LineWidth',2);
title('ԭʼ��ָ֤����ÿ�տ�������һ�����ͼ��','FontSize',12);
xlabel('����������(1990.12.19-2009.08.19)','FontSize',12);
ylabel('��һ����Ŀ�����','FontSize',12);
grid on;
% ��TS����ת��,�Է���libsvm����������ݸ�ʽҪ��
TS = TS';

% mapminmaxΪmatlab�Դ���ӳ�亯��
% ��tsx���й�һ��
[TSX,TSXps] = mapminmax(tsx,1,2);	
% ��TSX����ת��,�Է���libsvm����������ݸ�ʽҪ��
TSX = TSX';

%% ѡ��ع�Ԥ�������ѵ�SVM����c&g

% % ���Ƚ��д���ѡ��: 
% [bestmse,bestc,bestg] = SVMcgForRegress(TS,TSX,-8,8,-8,8);
% 
% % ��ӡ����ѡ����
% disp('��ӡ����ѡ����');
% str = sprintf( 'Best Cross Validation MSE = %g Best c = %g Best g = %g',bestmse,bestc,bestg);
% disp(str);
% 
% % ���ݴ���ѡ��Ľ��ͼ�ٽ��о�ϸѡ��: 
% [bestmse,bestc,bestg] = SVMcgForRegress(TS,TSX,-4,4,-4,4,3,0.5,0.5,0.05);
% 
% % ��ӡ��ϸѡ����
% disp('��ӡ��ϸѡ����');
% str = sprintf( 'Best Cross Validation MSE = %g Best c = %g Best g = %g',bestmse,bestc,bestg);
% disp(str);
% 
% %% ���ûع�Ԥ�������ѵĲ�������SVM����ѵ��
% cmd = ['-c ', num2str(bestc), ' -g ', num2str(bestg) , ' -s 3 -p 0.01'];

cmd = '-s 3 -t 2 -c 0.25 -g 2 -p 0.01';
model = libsvmtrain(TS,TSX,cmd);

%% SVM����ع�Ԥ��
[predict,mse,pro] = libsvmpredict(TS,TSX,model,'b -1');
% predict = mapminmax('reverse',predict',TSps);
% predict = predict';

% ��ӡ�ع���
str = sprintf( '������� MSE = %g ���ϵ�� R = %g%%',mse(2),mse(3)*100);
disp(str);

%% �������
figure;
plot(TS,'b-o');
hold on;
plot(predict,'r-^');
legend('ԭʼ����','�ع�Ԥ������');
hold off;
title('ԭʼ���ݺͻع�Ԥ�����ݶԱ�','FontSize',12);
xlabel('����������(1990.12.19-2009.08.19)','FontSize',12);
ylabel('������','FontSize',12);
grid on;

% figure;
% error = predict - TS';
% plot(error,'rd');
% title('���ͼ(predicted data - original data)','FontSize',12);
% xlabel('����������(1990.12.19-2009.08.19)','FontSize',12);
% ylabel('�����','FontSize',12);
% grid on;
% 
% figure;
% error = (predict - TS')./TS';
% plot(error,'rd');
% title('������ͼ(predicted data - original data)/original data','FontSize',12);
% xlabel('����������(1990.12.19-2009.08.19)','FontSize',12);
% ylabel('��������','FontSize',12);
% grid on;
snapnow;
toc;

%% �Ӻ��� SVMcgForRegress.m
function [mse,bestc,bestg] = SVMcgForRegress(train_label,train,cmin,cmax,gmin,gmax,v,cstep,gstep,msestep)
%SVMcg cross validation by faruto

%

% about the parameters of SVMcg 
if nargin < 10
    msestep = 0.06;
end
if nargin < 8
    cstep = 0.8;
    gstep = 0.8;
end
if nargin < 7
    v = 5;
end
if nargin < 5
    gmax = 8;
    gmin = -8;
end
if nargin < 3
    cmax = 8;
    cmin = -8;
end
% X:c Y:g cg:acc
[X,Y] = meshgrid(cmin:cstep:cmax,gmin:gstep:gmax);
[m,n] = size(X);
cg = zeros(m,n);

eps = 10^(-4);

bestc = 0;
bestg = 0;
mse = Inf;
basenum = 2;
for i = 1:m
    for j = 1:n
        cmd = ['-v ',num2str(v),' -c ',num2str( basenum^X(i,j) ),' -g ',num2str( basenum^Y(i,j) ),' -s 3 -p 0.1'];
        cg(i,j) = libsvmtrain(train_label, train, cmd);
        
        if cg(i,j) < mse
            mse = cg(i,j);
            bestc = basenum^X(i,j);
            bestg = basenum^Y(i,j);
        end
        
        if abs( cg(i,j)-mse )<=eps && bestc > basenum^X(i,j)
            mse = cg(i,j);
            bestc = basenum^X(i,j);
            bestg = basenum^Y(i,j);
        end
        
    end
end
% to draw the acc with different c & g
[cg,ps] = mapminmax(cg,0,1);
figure;
[C,h] = contour(X,Y,cg,0:msestep:0.5);
clabel(C,h,'FontSize',10,'Color','r');
xlabel('log2c','FontSize',12);
ylabel('log2g','FontSize',12);
firstline = 'SVR����ѡ����ͼ(�ȸ���ͼ)[GridSearchMethod]'; 
secondline = ['Best c=',num2str(bestc),' g=',num2str(bestg), ...
    ' CVmse=',num2str(mse)];
title({firstline;secondline},'Fontsize',12);
grid on;

figure;
meshc(X,Y,cg);
% mesh(X,Y,cg);
% surf(X,Y,cg);
axis([cmin,cmax,gmin,gmax,0,1]);
xlabel('log2c','FontSize',12);
ylabel('log2g','FontSize',12);
zlabel('MSE','FontSize',12);
firstline = 'SVR����ѡ����ͼ(3D��ͼ)[GridSearchMethod]'; 
secondline = ['Best c=',num2str(bestc),' g=',num2str(bestg), ...
    ' CVmse=',num2str(mse)];
title({firstline;secondline},'Fontsize',12);