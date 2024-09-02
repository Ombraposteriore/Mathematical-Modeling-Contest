%算法主程序
function GA4JSP205()
    clc;clear;
    %决定编码的主程序n*m，从1、2、、、n出现m次
    %dt=initDataFT06();%初始化数据
    dt=initDataFT10();
    [rows,cols]=size(dt);
    jobQty=rows;%作业数量
    machQty=cols/2;%机器数量
    pop=80;%群体规模
    chromes=createChromes(jobQty,machQty,pop);%初始化种群
    %种群第一代
    chrome1=chromes(1,:);
    sch=createSchedule(dt,chrome1);
    finishT=fitness(sch);
    
    %初始化算法参数
    crossRate=0.4; %交叉概率
    muteRate=0.6;%变异的比率
    maxGeneration=300;%循环代数
    nowGeneration=0;%当前代数
    optChrome=chrome1;%最优基因代数
    optFitness=finishT;%最优解
    
    %迭代循环程序-包含终止条件
    while nowGeneration<maxGeneration
        
        %遗传复制操作
        chromes=copyChromes(dt,chromes);
        %交叉操作
        %disp(['generation:' int2str(nowGeneration)]);
        %disp(chromes);
        chromes=crossChromes(chromes,crossRate);
        %变异操作
        chromes=muteChromesInv(chromes,muteRate);
        %判断是否替换最优解
        [nowFit,nowChrome]=findOptSolution(chromes,dt);
        if nowFit<optFitness
            optFitness=nowFit;
            optChrome=nowChrome;
        end    
        %optFitness
        nowGeneration=nowGeneration+1;
    end
    
    %显示最优结果
    disp(['最优结果：' int2str(optFitness)]);
    optSch=createSchedule(dt,optChrome)
    drawGant(optSch);%画出甘特图
end

%找出当前种群中最优解的染色体编码以及最优解的值
function [nowOptFit,nowOptChrome]=findOptSolution(inChromes,data)
    %step1：计算每条染色体的适应度值，放到数组中
    pop=size(inChromes,1);
    fitnesses=zeros(pop,1);
    for i=1:pop
        sch=createSchedule(data,inChromes(i,:));
        fit=fitness(sch);
        fitnesses(i)=fit;
    end
    %对适应度值进行排序
    [fit2,chromeId]=sortrows(fitnesses,1);
    nowOptFit=fit2(1); %适应度值中的第一个是最好的
    nowOptChrome=inChromes(chromeId(1),:);
end

%变异操作-逆序方法
function outChromes=muteChromesInv(chromes,muteRate)
    [pop,len]=size(chromes);
    for i=1:pop
        nowChrome=chromes(i,:);
        if rand()<muteRate
            %outPoint=randperm(len,1);
            %inPoint=randperm(len,1);
            points=sortrows(randperm(len,2)')';
            newChrome0=nowChrome(points(1):points(2)); %提取片段
            newChrome1=fliplr(newChrome0); %逆序操作
            chromes(i,points(1):points(2))=newChrome1;
        end  
        outChromes=chromes;
    end    
end

%变异操作-插入方法
function outChromes=muteChromesInster(chromes,muteRate)
    [pop,len]=size(chromes);
    for i=1:pop
        nowChrome=chromes(i,:); %拿出当前行
        if rand()<muteRate
            %outPoint=randperm(len,1);
            %inPoint=randperm(len,1);
            points=randperm(len,2); %随机取两点，points(1)为插入点，points(2)为被插入点
            if points(1)<points(2)
                insertChrome=nowChrome(points(1):points(2)-1);
                newChrome=circshift(insertChrome,-1);
                chromes(i,points(1):points(2)-1)=newChrome;
            else
                insertChrome=nowChrome(points(2)+1:points(1));
                newChrome=circshift(insertChrome,1);
                chromes(i,points(2)+1:points(1))=newChrome;
            end
            chromes(i,:);
        end  
        outChromes=chromes;
    end    
end

 %部分映射交叉操作
function outChromes=crossChromes(chromes,crossRate)
    pop=size(chromes,1);
    cols=size(chromes,2);
    for i=1:2:pop
        if rand()<crossRate %进行交叉操作
            parent1=chromes(i,:);
            parent2=chromes(i+1,:);
            p=sortrows(randperm(cols,2)')'; %随机生成两个数后按列排序
            %p=min(p):max(p);
            %两个交换片段
            crossP1=parent1(p(1):p(2));
            crossP2=parent2(p(1):p(2));
            son1=parent1;
            son1(p(1):p(2))=crossP2;
            son2=parent2;
            son2(p(1):p(2))=crossP1;
            %子片段对比
            crossLen=size(crossP1,2);
            for j=crossLen:-1:1
                midCode=crossP1(j);
                for k=1:size(crossP2,2)
                    if crossP2(k)==midCode
                        crossP1(j)=[];
                        crossP2(k)=[];
                        break;
                    end
                end
            end 
            %染色体编码置换
            repeatNum=size(crossP1,2); %重复编码片段
            if repeatNum>0
                %对子代1的有效性置换
                for j=1:repeatNum
                    midCode=crossP2(j);
                    if p(1)==1 %判断片段后是否重复
                        for k=p(2)+1:cols
                            if son1(k)==midCode
                                son1(k)=crossP1(j);
                                break;
                            end
                        end
                    else
                        if p(2)==cols
                           for k=1:p(1)-1
                                if son1(k)==midCode
                                    son1(k)=crossP1(j);
                                    break;
                                end
                            end 
                        else
                            getIt=0;
                            for k=1:p(1)-1
                                if son1(k)==midCode
                                    son1(k)=crossP1(j);
                                    getIt=1; %找到了为1
                                    break;
                                end
                            end
                            if getIt==0
                                for k=p(2)+1:cols
                                    if son1(k)==midCode
                                        son1(k)=crossP1(j);
                                        break;
                                    end
                                 end
                            end    
                        end    
                    end
                end    
                %对子代2的有效性置换 
                for j=1:repeatNum
                    midCode=crossP1(j);
                    if p(1)==1
                        for k=p(2)+1:cols
                            if son2(k)==midCode
                                son2(k)=crossP2(j);
                                break;
                            end
                        end
                    else
                        if p(2)==cols
                           for k=1:p(1)-1
                                if son2(k)==midCode
                                    son2(k)=crossP2(j);
                                    break;
                                end
                            end 
                        else
                            getIt=0;
                            for k=1:p(1)-1
                                if son2(k)==midCode
                                    son2(k)=crossP2(j);
                                    getIt=1;
                                    break;
                                end
                            end
                            if getIt==0
                                for k=p(2)+1:cols
                                    if son2(k)==midCode
                                        son2(k)=crossP2(j);
                                        break;
                                    end
                                 end
                            end    
                        end    
                    end
                end  
            end 
            chromes(i,:)=son1;
            chromes(i+1,:)=son2;
        end    
    end    
    outChromes=chromes;
end

%轮盘赌策略的复制操作
function outChromes=copyChromes(data,inChromes)
    outChromes=zeros(size(inChromes));
    %1：计算每条染色体的适应度值，放到数组中
    pop=size(inChromes,1);%输入染色体的行数
    fitnesses=zeros(pop,1); %存储每条染色体的适应度值
    for i=1:pop 
        sch=createSchedule(data,inChromes(i,:));%计算种群当中的第i条染色体的调度方案
        fit=fitness(sch); %计算最优值
        fitnesses(i)=fit;
    end
    [fit2,chromeId]=sortrows(fitnesses,1); %fit2为排序后的数值序列，chromeId为行号的排序序列
    maxVal=max(fit2)*1.2; %定义一个参考值
    fit3=maxVal-fit2; %排序后的适应度参考值
    sumFit=sum(fit3); %求和
    oneFit=fit3/sumFit; %归一化处理
    finalFit=zeros(pop,1);
    for i=1:pop %累加概率值，当前值为前面数值之和
        finalFit(i)=sum(oneFit(1:i));
    end
    %从现有pop个染色体中轮盘赌选择新成新的染色体
    for i=1:pop
        midVal=rand();
        if midVal<finalFit(1) 
            outChromes(i,:)=inChromes(chromeId(1),:);
        else
            for j=1:pop-1
                if midVal>=finalFit(j)&&midVal<finalFit(j+1)
                    outChromes(i,:)=inChromes(chromeId(j+1),:);
                    break;
                end  
            end 
        end    
    end 
    
end

%获取调度方案的最大完工时间
function finishTime=fitness(schedule)
    finishTime=max(schedule(:,5));
end

%根据调度方案绘制甘特图
function drawGant(schedule)
    rows=size(schedule,1);
    maxMachId=max(schedule(:,2));
    jobQty=max(schedule(:,1));
    mycolor=rand(jobQty,3);%随机生成颜色
    figure;
    ylim([0 maxMachId+1]); %设置y轴上下限
    title('甘特图');
    xlabel('时间');
    ylabel('机器');
    for i=1:rows
        x=schedule(i,4:5); %开始时间到结束时间
        y=[schedule(i,2) schedule(i,2)]; %机器Id
        line(x,y,'lineWidth',16,'color',mycolor(schedule(i,1),:)); 
        procId=schedule(i,3);
        jobId=schedule(i,1);
        txt=['[' int2str(jobId) '-' int2str(procId) ']']; %当前工序编号和作业编号
        text(mean(x)-1,y(1),txt);
    end
end

%将编码转换成具体的调度方案
function schedule=createSchedule(data,chrome)
    jobQty=size(data,1);  %返回数组第一个维度的大小，行数
    machQty=size(data,2)/2; %返回数组的列数,两列一个设备，故除2
    %用这个获得行数列数也行： [rows,cols]=size(data)
    schedule=zeros(jobQty*machQty,5);
    %定义中间数组
    %当前作业的开始时间
    jobCanStartTime=zeros(1,jobQty);
    %当前拟排产的工序号
    jobProcessId=ones(1,jobQty);
    %获取编码中当前作业号，i位置的数字
    for i=1:jobQty*machQty
        nowJobId=chrome(i);
        nowProcId=jobProcessId(nowJobId);
        nowMachId=data(nowJobId,2*nowProcId-1);
        nowProcTime=data(nowJobId,2*nowProcId);
        machSch=schedule(schedule(:,2)==nowMachId,:);%取出第二列等于当前机器Id的行
        jobCanST=jobCanStartTime(nowJobId);
        if size(machSch,1)==0 %设备还没有安排作业
            startTime=jobCanStartTime(nowJobId);
            endTime=startTime+nowProcTime;%结束时间为开始时间+工序耗时
        else %设备已经安排了作业
            machSch=sortrows(machSch,4);%按第四列开始时间进行行排序
            rows=size(machSch,1);
            %处理第一行已排作业，检查是否能够将当前作业排到该作业之前
            done=0;
            if jobCanST<machSch(1,4) %第一行开始排
                if machSch(1,4)-jobCanST>=nowProcTime %工时-可开工时间大于工序耗时，则可以排到前面
                    startTime=jobCanST;
                    endTime=startTime+nowProcTime;
                    done=1;
                end    
            end    
            if done==0
                for j=2:rows %从第二行开始循环判断，是否能够排到其他作业之前
                    if jobCanST<machSch(j,4)
                        if machSch(j,4)-max(jobCanST,machSch(j-1,5))>=nowProcTime
                            startTime=max(jobCanST,machSch(j-1,5));
                            endTime=startTime+nowProcTime;
                            done=1;
                            break;
                        end
                    end
                end   
            end 
            if done==0 %表示该作业不能排到已有作业之前，则排到最后
                startTime=max(jobCanST,machSch(rows,5));
                endTime=startTime+nowProcTime;
            end 
        end
        schedule(i,1)=nowJobId;%第一列作业Id
        schedule(i,2)=nowMachId;%第二列机器Id
        schedule(i,3)=nowProcId;%第三列工序Id
        schedule(i,4)=startTime;%第四列开始时间
        schedule(i,5)=endTime;%第五列结束时间
        jobCanStartTime(nowJobId)=endTime; %下一个开始的时间等于上一个的结束时间
        jobProcessId(nowJobId)=jobProcessId(nowJobId)+1;
    end    
end

%初始化工时设备数据
function data = initDataFT06()
%FT06的数据，一行表示一个作业，共六个作业、六个机器数（机器号下表从0开始）
%奇数列1、3、5、7、9、11表示该工序使用的机器号
%偶数列2、4、6、8、10、12表示在该机器上加工所用时间
    data=[ 2  1  0  3  1  6  3  7  5  3  4  6;
 1  8  2  5  4 10  5 10  0 10  3  4;
 2  5  3  4  5  8  0  9  1  1  4  7;
 1  5  0  5  2  5  3  3  4  8  5  9;
 2  9  1  3  4  5  5  4  0  3  3  1;
 1  3  3  3  5  9  0 10  4  4  2  1];
    fullCols=size(data,2); %获取列数
    cols=1:2:fullCols;
    %将data奇数列设备编码由原来的0开始，转换为从1开始
    data(:,cols)=data(:,cols)+1;
end

function data = initDataFT10()
%奇数列1、3、5、7、9、11表示该工序使用的机器号
%偶数列2、4、6、8、10、12表示在该机器上加工所用时间
    data=[ 0 29 1 78 2  9 3 36 4 49 5 11 6 62 7 56 8 44 9 21;
 0 43 2 90 4 75 9 11 3 69 1 28 6 46 5 46 7 72 8 30;
 1 91 0 85 3 39 2 74 8 90 5 10 7 12 6 89 9 45 4 33;
 1 81 2 95 0 71 4 99 6  9 8 52 7 85 3 98 9 22 5 43;
 2 14 0  6 1 22 5 61 3 26 4 69 8 21 7 49 9 72 6 53;
 2 84 1  2 5 52 3 95 8 48 9 72 0 47 6 65 4  6 7 25;
 1 46 0 37 3 61 2 13 6 32 5 21 9 32 8 89 7 30 4 55;
 2 31 0 86 1 46 5 74 4 32 6 88 8 19 9 48 7 36 3 79;
 0 76 1 69 3 76 5 51 2 85 9 11 6 40 7 89 4 26 8 74;
 1 85 0 13 2 61 6  7 8 64 9 76 5 47 3 52 4 90 7 45];
    fullCols=size(data,2);
    cols=1:2:fullCols;
    %将data奇数列设备编码由原来的0开始，转换为从1开始
    data(:,cols)=data(:,cols)+1;
end

%根据种群数量和作业数量和设备数量生成初始种群
function chromes=createChromes(jobQty,machQty,pop)
    chromes=zeros(pop,jobQty*machQty);
    for i=1:pop
        %每一行chromes都由createChrome生成
        chromes(i,:)=createChrome(jobQty,machQty);
    end
end

%生成编码：传两个参数（1）jobQty，（2）machQty
%jobQty:作业数量  machQty:机器数量
function chrome=createChrome(jobQty,machQty)
    a=1:jobQty;
    chrome=a;
    for i=1:machQty
        chrome=[chrome,a];
    end
    %将编码乱序排列-生成作业码的乱序
    b=randperm(jobQty*machQty);
    chrome=chrome(b);
end