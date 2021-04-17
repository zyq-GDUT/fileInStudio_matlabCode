%%%%%%%%%%%%%%%%%免疫算法求函数极值%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%初始化%%%%%%%%%%%%%%%%%%%%%%%%%
clear all;                                %清除所有变量
close all;                                %清图
%clc;                                      %清屏
global fes
MAXRUN = 20;                %独立测试次数
D=10;                                     %免疫个体维数
NP=100;                                   %免疫个体数目
Xs=20;                                    %取值上限
Xx=-20;                                   %取值下限
G=1000;                                    %最大免疫代数
pm=0.7;                                   %变异概率
alfa=1;                                   %激励度系数
belta=1;                                  %激励度系数
detas=0.2;                                %相似度阈值
Ncl=10;                                   %克隆个数
deta0=1*Xs;                               %邻域范围初值
error = 10^-6;                  %获得解的误差值

tallrecord = [];
bestfrecord = [];
meantrace = [];
bestGenRecord = [];
bestFEsRecord = [];
bestTimeRecord = [];
okNum = 0;
allfile = fopen('total.txt','a');
for run = 1:MAXRUN
    fes = 0;
    %%%%%%%%%%%%%%%%%%%%%%%初始种群%%%%%%%%%%%%%%%%%%%%%%%%
    f=rand(D,NP)*(Xs-Xx)+Xx;            %NP个列向量个体
    for np=1:NP
        FIT(np)=func1(f(:,np));         %亲和度 = 目标函数值
    end
    gbestFIT = FIT(1);      %随机选第一个个体的目标函数值充当初始的至今最优解值 +++++++++
    gbextX = f(:,1);        %随机选第一个个体作为至今最优解+++++++++++
    bestGen = 0;
    bestFEs = fes;
    besttime = 0;
    %%%%%%%%%%%%%%%%%计算个体浓度和激励度%%%%%%%%%%%%%%%%%%%
    for np=1:NP
        for j=1:NP
            nd(j)=sqrt(sum((f(:,np)-f(:,j)).^2)); %抗体间亲和度P63（4.3）
            if nd(j)<detas          %P62（4.2）中的相似度阈值（=0.2）
                nd(j)=1;            %个体i和个体j相似
            else
                nd(j)=0;            %个体i和个体j距离比较远，不近似
            end
        end
        ND(np)=sum(nd)/NP;          %抗体浓度P62（4.1），越相似浓度越大
    end
    FIT = alfa*FIT- belta*ND; %基于目标函数值和抗体浓度得到个体激励度P63（4.6）
    %%%%%%%%%%%%%%%%%%%激励度按升序排列%%%%%%%%%%%%%%%%%%%%%%
    [SortFIT,Index]=sort(FIT);
    Sortf=f(:,Index);               %本例值越小越好
    
    onerunfile = fopen(['F1_run' num2str(run) '.txt'],'w'); %w 改为a表示添加   
    t0 = clock;
    tall = 0;                       %单次测试总用时
    gen=0;                          %免疫代数
    flagOK = false;
    %%%%%%%%%%%%%%%%%%%%%%%%免疫循环%%%%%%%%%%%%%%%%%%%%%%%%
    while gen<G
        for i=1:NP/2    %选激励度前NP/2个体进行免疫操作%%%%%%%%%%%
            a=Sortf(:,i);         %第i个个体的编码a
            Na=repmat(a,1,Ncl);   %克隆编码Nc1份，即1行Nc1列（Nc1=10）
            deta=deta0/gen;       %变异邻域范围值随代数gen变小
            for j=1:Ncl           %对克隆的每一个副本
                for ii=1:D
                    %%%%%%%%%%%%%%%%%变异%%%%%%%%%%%%%%%%%%%
                    if rand<pm    %如果满足变异概率
                        Na(ii,j)=Na(ii,j)+(rand-0.5)*deta; %邻域内的局部搜索P64（4.9）
                    end
                    %%%%%%%%%%%%%%边界条件处理%%%%%%%%%%%%%%%
                    if (Na(ii,j)>Xs)  ||  (Na(ii,j)<Xx)
                        Na(ii,j)=rand * (Xs-Xx)+Xx;     %随机值替换
                    end
                end
            end
            
            Na(:,1)=Sortf(:,i); %保留克隆源个体，此处做法：原个体取代克隆的第一个变异副本
            
            %%%%%%%%%%克隆抑制，保留亲和度最高的个体%%%%%%%%%%
            for j=1:Ncl
                NaFIT(j)=func1(Na(:,j));    %对每个克隆副本（第一个是原个体）函数评价
            end
            [NaSortFIT,Index]=sort(NaFIT);
            aFIT(i)=NaSortFIT(1);           %最优的函数评价值保存在aFIT(i)
            NaSortf=Na(:,Index);
            af(:,i)=NaSortf(:,1);           %函数评价最优个体保存在af的第i列
            
            if aFIT(i) < gbestFIT             %比较是否找到更优的至今最优解+++++
                gbestFIT = aFIT(i);           %更新至今最优解值
                gbextX = af(:,i);          %更新至今最优解
                bestGen = gen;                 %首次获得至今最优解值对应的平均代数
                bestFEs = fes;                 %首次获得至今最优解值对应的函数评价次数
                besttime = tall + etime(clock,t0);  %对应的用时
                if(~flagOK && gbestFIT < error) %首次找到满足精度error内的解
                    flagOK = true;
                    okNum = okNum +1;
                end
            end
        end
        %%%%%%%%%%%%%%%%%%%%免疫种群激励度%%%%%%%%%%%%%%%%%%%
        for np=1:NP/2               %一半种群
            for j=1:NP/2
                nda(j)=sqrt(sum((af(:,np)-af(:,j)).^2)); %抗体间亲和度P63（4.3）
                if nda(j)<detas      %P62（4.2）中的相似度阈值（=0.2）
                    nda(j)=1;
                else
                    nda(j)=0;
                end
            end
            aND(np)=sum(nda)/NP/2;  %此处为一半种群
        end
        aFIT =  alfa*aFIT-  belta*aND;  %得到一半种群的抗体激励度
        %%%%%%%%%%%%%%%%%%%%%%%种群刷新%%%%%%%%%%%%%%%%%%%%%%%
        bf=rand(D,NP/2)*(Xs-Xx)+Xx;     %另一半种群随机生成
        for np=1:NP/2
            bFIT(np)=func1(bf(:,np));
        end
        
        [ibestFIT,irank] = min(bFIT);      %本次迭代中函数值最优的值+++++++++
        if ibestFIT < gbestFIT             %比较是否找到更优的至今最优解+++++
            gbestFIT = ibestFIT;           %更新至今最优解值
            gbextX = bf(:,irank);          %更新至今最优解
            bestGen = gen;                 %首次获得至今最优解值对应的平均代数
            bestFEs = fes;                 %首次获得至今最优解值对应的函数评价次数
            besttime = tall + etime(clock,t0);  %对应的用时
            if(~flagOK && gbestFIT < error) %首次找到满足精度error内的解
                flagOK = true;
                okNum = okNum +1;
            end
        end
        
        %%%%%%%%%%%%%%%%%%%新生成种群激励度%%%%%%%%%%%%%%%%%%%%
        for np=1:NP/2               %一半种群
            for j=1:NP/2
                ndc(j)=sqrt(sum((bf(:,np)-bf(:,j)).^2));
                if ndc(j)<detas
                    ndc(j)=1;
                else
                    ndc(j)=0;
                end
            end
            bND(np)=sum(ndc)/NP/2;    %此处为一半种群
        end
        bFIT =  alfa*bFIT-  belta*bND; %得到另一半种群的抗体激励度
        %%%%%%%%%%%%%%免疫种群与新生种群合并%%%%%%%%%%%%%%%%%%%
        f1=[af,bf];                   %两半种群个体列合并
        FIT1=[aFIT,bFIT];             %两半种群抗体激励度合并
        [SortFIT,Index]=sort(FIT1);   %整个种群按抗体激励度排序
        Sortf=f1(:,Index);            %获得排序后的新种群
        gen=gen+1;
        trace(gen)= gbestFIT; %func1(Sortf(:,1)); %记录最优激励度抗体对应的函数评价值++++
        
        tall = tall + etime(clock,t0);
        sprintf('%d\t%d\t%g\r\n',gen,fes,trace(gen)); %%%%%%out%%%%%%%%%%%
        fprintf(onerunfile,'%d\t%d\t%g\r\n',gen,fes,trace(gen));
        t0 = clock;
    end
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
end

meantrace = mean(meantrace);     %按列取均值，获得平均收敛曲线
meanfile = fopen('F1runMean.txt','w');
for i=1:gen
    fprintf(meanfile,'%g\r\n',meantrace(i));
end
fclose(meanfile);

bestfrecord;
fprintf(allfile,'%s\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\r\n','F1',min(bestfrecord),max(bestfrecord),...
    mean(bestfrecord),std(bestfrecord),...
    mean(bestGenRecord),mean(bestFEsRecord),mean(bestTimeRecord),...
    mean(tallrecord),okNum);
fclose(allfile);

%%%%%%%%%%%%%%%%%%%%%%%输出优化结果%%%%%%%%%%%%%%%%%%%%%%%%
% Bestf=Sortf(:,1);          %最优个体变量
% trace(end);                %最优值，end代表取trace向量最后一个数
% figure,plot(trace)
% xlabel('迭代次数')
% ylabel('目标函数值')
% title('亲和度进化曲线')