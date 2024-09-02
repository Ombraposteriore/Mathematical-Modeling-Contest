
clear
clc
close all
tic
%% ��importdata�����������ȡ�ļ�
c101=importdata('c101.txt');
cap=500;                                                        %�������װ����
%% ��ȡ������Ϣ
E=c101(1,5);                                                    %��������ʱ�䴰��ʼʱ��
L=c101(1,6);                                                    %��������ʱ�䴰����ʱ��
vertexs=c101(:,2:3);                                            %���е������x��y
customer=vertexs(2:end,:);                                       %�˿�����
cusnum=size(customer,1);                                         %�˿���
v_num=2;                                                        %�������ʹ����Ŀ
demands=c101(2:end,4);                                          %������
a=c101(2:end,5);                                                %�˿�ʱ�䴰��ʼʱ��[a[i],b[i]]
b=c101(2:end,6);                                                %�˿�ʱ�䴰����ʱ��[a[i],b[i]]
s=c101(2:end,7);                                                %�ͻ���ķ���ʱ��
% ������еĽڵ�
nodes = 16;
% ������еı�
edges=[1    1	0	0	1	0	1	0	0	0	1	0	0	0	0	0;
1	1	1	0	1	0	0	0	0	0	0	0	0	0	0	0;
0	1	1	1	1	1	0	0	0	0	0	0	0	0	0	0;
0	0	1	1	0	1	0	0	0	1	0	0	0	0	0	0;
1	1	1	0	1	1	1	1	1	0	0	0	0	0	0	0;
0	0	1	1	1	1	0	0	1	1	0	0	0	0	0	0;
1	0	0	0	1	0	1	1	0	0	1	1	0	0	0	0;
0	0	0	0	1	0	1	1	1	0	0	1	1	0	1	0;
0	0	0	0	1	1	0	1	1	1	0	1	1	1	1	1;
0	0	0	1	0	1	0	0	1	1	0	0	0	1	0	1;
1	0	0	0	0	0	1	0	0	0	1	1	0	0	0	0;
0	0	0	0	0	0	1	1	1	0	1	1	1	0	1	0;
0	0	0	0	0	0	0	1	1	0	0	1	1	1	1	1;
0	0	0	0	0	0	0	0	1	1	0	0	1	1	1	1;
0	0	0	0	0	0	0	1	1	0	0	1	1	1	1	0;
0	0	0	0	0	0	0	0	1	1	0	0	1	1	0	1
];
% �ͻ�λ�ú����������� ��Ĭ�ϵ�һ��Ϊ�ֿ�
locations =[20.20818722,46.55925669;30.86325116,42.42106760;43.61357686,38.9299278;
    35.00118717,36.98696184;26.99611373,17.4952427;46.18764943,42.83843589;
    40.6920241,31.94898546;32.50662006,10.62166377;45.88614135,13.7874543;
    24.72317638,30.59033166;23.50924028,48.09930553;34.75421424,43.72244785;
    13.4758176,47.40048984;25.01990127,48.9365234;13.4758176,47.40048984;25.01990127,48.9365234];
%   ����������
dist = zeros(nodes, nodes);
for i = 1:nodes
    for j = 1:nodes
        if edges(i,j)==1  % �ж��Ƿ���ֱ������
            dist(i,j) = norm(locations(i,:) - locations(j,:));
        else
            dist(i,j) = 100000;  % �������ֱ�����ӣ���������Ϊ���޴�
        end
    end
end 
% h=pdist(vertexs);
% dist=squareform(h);                                             %��������������ǹ�ϵ�����þ����ʾ����c[i][j]=dist[i][j]
%% �Ŵ��㷨��������
alpha=1000;                                                       %Υ��������Լ���ĳͷ�����ϵ��
belta=0;                                                      %Υ��ʱ�䴰Լ���ĳͷ�����ϵ��
NIND=300;                                                       %��Ⱥ��С
MAXGEN=800;                                                     %��������
Pc=0.91;                                                         %�������
Pm=0.01;                                                        %�������
GGAP=0.9;                                                       %����(Generation gap)
N=cusnum+v_num-1;                                               %Ⱦɫ�峤��=�˿���Ŀ+�������ʹ����Ŀ-1
%% ��ʼ����Ⱥ
init_vc=init(cusnum,a,demands,cap);                             %�����ʼ��
Chrom=InitPopCW(NIND,N,cusnum,init_vc);
%% ���������·�ߺ��ܾ���
disp('��ʼ��Ⱥ�е�һ�����ֵ:')
[VC,NV,TD,violate_num,violate_cus]=decode(Chrom(1,:),cusnum,cap,demands,a,b,L,s,dist);
% disp(['�ܾ��룺',num2str(TD)]);
disp(['����ʹ����Ŀ��',num2str(NV),'��������ʻ�ܾ��룺',num2str(TD),'��Υ��Լ��·����Ŀ��',num2str(violate_num),'��Υ��Լ���˿���Ŀ��',num2str(violate_cus)]);
disp('~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~')
%% �Ż�
gen=1;
figure;
hold on;box on
xlim([0,MAXGEN])
title('�Ż�����')
xlabel('����')
ylabel('����ֵ')
ObjV=calObj(Chrom,cusnum,cap,demands,a,b,L,s,dist,alpha,belta);             %������ȺĿ�꺯��ֵ
preObjV=min(ObjV);
while gen<=MAXGEN
    %% ������Ӧ��
    ObjV=calObj(Chrom,cusnum,cap,demands,a,b,L,s,dist,alpha,belta);             %������ȺĿ�꺯��ֵ
    line([gen-1,gen],[preObjV,min(ObjV)]);pause(0.0001)
    preObjV=min(ObjV);
    FitnV=Fitness(ObjV);
    %% ѡ��
    SelCh=Select(Chrom,FitnV,GGAP);
    %% OX�������
    SelCh=Recombin(SelCh,Pc);
    %% ����
    SelCh=Mutate(SelCh,Pm);
    %% �ֲ���������
    SelCh=LocalSearch(SelCh,cusnum,cap,demands,a,b,L,s,dist,alpha,belta);
    %% �ز����Ӵ�������Ⱥ
    Chrom=Reins(Chrom,SelCh,ObjV);
    %% ɾ����Ⱥ���ظ����壬������ɾ���ĸ���
    Chrom=deal_Repeat(Chrom);
    %% ��ӡ��ǰ���Ž�
    ObjV=calObj(Chrom,cusnum,cap,demands,a,b,L,s,dist,alpha,belta);             %������ȺĿ�꺯��ֵ
    [minObjV,minInd]=min(ObjV);
    disp(['��',num2str(gen),'�����Ž�:'])
    [bestVC,bestNV,bestTD,best_vionum,best_viocus]=decode(Chrom(minInd(1),:),cusnum,cap,demands,a,b,L,s,dist);
    disp(['����ʹ����Ŀ��',num2str(bestNV),'��������ʻ�ܾ��룺',num2str(bestTD),'��Υ��Լ��·����Ŀ��',num2str(best_vionum),'��Υ��Լ���˿���Ŀ��',num2str(best_viocus)]);
    fprintf('\n')
    %% ���µ�������
    gen=gen+1 ;
end
%% �������Ž��·��ͼ
ObjV=calObj(Chrom,cusnum,cap,demands,a,b,L,s,dist,alpha,belta);             %������ȺĿ�꺯��ֵ
[minObjV,minInd]=min(ObjV);
%% ������Ž��·�ߺ��ܾ���
disp('���Ž�:')
bestChrom=Chrom(minInd(1),:);
[bestVC,bestNV,bestTD,best_vionum,best_viocus]=decode(bestChrom,cusnum,cap,demands,a,b,L,s,dist);
disp(['����ʹ����Ŀ��',num2str(bestNV),'��������ʻ�ܾ��룺',num2str(bestTD),'��Υ��Լ��·����Ŀ��',num2str(best_vionum),'��Υ��Լ���˿���Ŀ��',num2str(best_viocus)]);
disp('-------------------------------------------------------------')
%% �ж����Ž��Ƿ�����ʱ�䴰Լ����������Լ����0��ʾΥ��Լ����1��ʾ����ȫ��Լ��
flag=Judge(bestVC,cap,demands,a,b,L,s,dist);
%% ������Ž����Ƿ����Ԫ�ض�ʧ���������ʧԪ�أ����û����Ϊ��
DEL=Judge_Del(bestVC);
%% ��������·��ͼ
draw_Best(bestVC,vertexs);
save c101.mat
toc