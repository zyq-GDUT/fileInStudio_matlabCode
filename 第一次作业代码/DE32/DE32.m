%% 初始化 %%
clear all;
close all;
clc;
NP=50;
D=10;
G=2000;
F=0.5;
CR=0.1;
Xs=20;
Xx=-20;
global x trace Ob ;

%% 输出数据的定义 %%
error=10^-6;
global fes;
MAXRUN =20;
tallrecord = [];
bestfrecord = [];
meantrace = [];
bestGenRecord = [];
bestFEsRecord = [];
bestTimeRecord = [];
okNum = 0;
allfile = fopen('total.txt','a');
timecount=zeros(1,MAXRUN);
%% 开始测试 %%
for run =1:MAXRUN
    tt1=clock;
    fes = 0;
    t0 = clock;
    %% 初始化种群 %%
    x=zeros(D,NP);
    v=zeros(D,NP);
    u=zeros(D,NP);
    x=rand(D,NP)*(Xs-Xx)+Xx;
    %% 计算适应度 %%
    for m =1:NP
        Ob(m)=funcl(x(:,m));
        fes=fes+1;
    end
    trace(1)=min(Ob);
    %% 输出数据准备 %%
    gbestFIT = Ob(1);  %随机选第一个个体的目标函数值充当初始的至今最优解值
    bestGen = 0;
    bestFEs = fes;
    besttime = 0;
    onerunfile = fopen(['F1_run' num2str(run) '.txt'],'w'); %w 改为a表示添加   
    tall = 0;        %单次测试总用时
    gen=0;           %免疫代数
    flagOK = false;  %到达误差内的最优值时转为true
    %% 差分循环 %%
    for gen=1:G
        v = mutation32(x,NP,F);
        u = crosser32(D,CR,v,x,NP,Xx,Xs);
        selection32(NP,u,gen)
        if trace(gen+1) < gbestFIT             %比较是否找到更优的至今最优解+++++
            gbestFIT = trace(gen+1);           %更新至今最优解值
            bestGen = gen+1;                   %首次获得至今最优解值对应的平均代数
            bestFEs = fes;                     %首次获得至今最优解值对应的函数评价次数
            besttime = tall + etime(clock,t0); %对应的用时
            if(~flagOK && gbestFIT < error)    %首次找到满足精度error内的解
                flagOK = true;
                okNum = okNum +1;
            end
        end
        tall = tall + etime(clock,t0);
        sprintf('%d\t%d\t%g\r\n',gen,fes,trace(gen)); %%%%%%out%%%%%%%%%%%
        fprintf(onerunfile,'%d\t%d\t%g\r\n',gen,fes,trace(gen));
        t0 = clock;
    end
    tt2=clock;
    fclose(onerunfile);
    tallrecord = [tallrecord,tall];
    bestfrecord = [bestfrecord, trace(gen)];
    funfile = fopen('F1total.txt','a'); %w 改为a表示添加
    fprintf('%d\t%g\t%g\r\n',run,trace(gen),tall)
    fprintf(funfile,'%d\t%g\t%g\t%g\t%g\t%g\r\n',run,trace(gen),bestGen,bestFEs,besttime,tall);
    fclose(funfile);
    meantrace = [meantrace; trace];  %按行合并迭代收敛行记录
    bestGenRecord = [bestGenRecord,bestGen];
    bestFEsRecord = [bestFEsRecord,bestFEs];
    bestTimeRecord = [bestTimeRecord,besttime];
    timecount(run)=etime(tt2,tt1);
end
 timecount
meantrace = mean(meantrace);     %按列取均值，获得平均收敛曲线
meanfile = fopen('F1runMean.txt','w');
for i=1:gen
    fprintf(meanfile,'%g\r\n',meantrace(i));
end
fclose(meanfile);

fprintf(allfile,'%s\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\r\n','F1',min(bestfrecord),max(bestfrecord),...
    mean(bestfrecord),std(bestfrecord),...
    mean(bestGenRecord),mean(bestFEsRecord),mean(bestTimeRecord),...
    mean(tallrecord),okNum);
fclose(allfile);
%% 画图 %%
%figure
%plot(trace);
%xlabel('迭代次数')
%ylabel('目标函数值')
%title('DE目标函数曲线')

