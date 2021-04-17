clear all;
close all;
clc;
%%初始化参数
NP=100;
Pc=0.1;
Pm=0.3;
evaluateCount=10^5;%最大迭代次数
G=ceil(evaluateCount/NP);%迭代次数
wonderful=15609.5;
error=10^-6;
MAXRUNS=20;
%各个城市的坐标，用于计算相互之间的距离%
C=[1304 2312;3639 1315;4177 2244;3712 1399;3488 1535;3326 1556;...
3238 1229;4196 1044;4312  790;4386  570;3007 1970;2562 1756;...
2788 1491;2381 1676;1332  695;3715 1678;3918 2179;4061 2370;...
3780 2212;3676 2578;4029 2838;4263 2931;3429 1908;3507 2376;...
3394 2643;3439 3201;2935 3240;3140 3550;2545 2357;2778 2826;...
2370 2975];                
city_num=size(C,1);
Dist=getDistance(C,city_num);%获取任意两个城市的距离%
run=0;
avgBestPerGen=zeros(G,1);%用来记录每次测试的平均每代最优
singleTestTime=zeros(MAXRUNS,1);%用来记录每次测试的时间
getBestTime=zeros(MAXRUNS,1);%用来记录每次测试中获取到最优值的时间
BestFitPerTest=zeros(MAXRUNS,1);%用来记录每次测试的最优适应度
genNum=zeros(MAXRUNS,1);%用来记录每次测试的最优代数
while run<MAXRUNS
    time1=clock;
    gen=1;
    ec=0;
    pop=initialPop(NP,city_num);
    bestTimePerGen=zeros(G,1);%每代获取到最优的时间
    Fit_best_perGen=zeros(G,1);
    while ec<evaluateCount    
        fit=fitnessEvaluation(pop,NP,Dist);
        [bestfit,index]=min(fit);
        bestPop=pop(index(1,1),:);
        getbestTime=clock;%%
        bestTimePerGen(gen,1)=etime(getbestTime,time1);%%
        Fit_best_perGen(gen,1)=bestfit;%% 
        ec=ec+NP;
        gen=gen+1;
        pop=roulette(pop,fit,NP,city_num);
        pop=cross(pop,NP,Pc,city_num,bestPop);
        pop=mutate(pop,NP,Pm,city_num);
        pop(1,:)=bestPop;
    end
    time2=clock;
    singleTestTime(run+1)=etime(time2,time1);
    [best,index]=min(Fit_best_perGen);
    getBestTime(run+1,1)=bestTimePerGen(index,1);
    BestFitPerTest(run+1,1)=best;
    genNum(run+1,1)=index;
    avgBestPerGen=avgBestPerGen+Fit_best_perGen;
    
    run=run+1;
end
avgBestPerGen= avgBestPerGen/MAXRUNS;
meanfit=mean(BestFitPerTest);
std_bestfit=std(BestFitPerTest);
bestfit=min(BestFitPerTest);%所有测试中最优值的最优值
worstfit=max(BestFitPerTest);%所有测试中最优值中的最差值
meangen=mean(genNum);%每次测试首次获得最优适应度的代数
meanruntime=mean(singleTestTime);%单次测试的平均时间
avg_evaluate_num=meangen*NP;
avg_besttime=mean(getBestTime);
success_count=0;
for i=1:MAXRUNS
    if (BestFitPerTest(i)-wonderful)<=error
        success_count=success_count+1;
    end
end
avg_success_count=success_count/MAXRUNS;
file=fopen('testdata.txt','w');
fprintf(file,'算法名：遗传算法解决TSP问题\n\n测试次数：%d\t迭代次数：%d\n\n种群规模：%d\t交叉概率：%.3f\t变异概率:%.3f\n\n',MAXRUNS,G,NP,Pc,Pm);
fprintf(file,'最优适应度的标准差：%.9g\n\n',std_bestfit);
fprintf(file,'最优适应度的最优值：%.9g\t',bestfit);
fprintf(file,'最优适应度的最差值：%.9g\n\n',worstfit);
fprintf(file,'获得最优适应度的平均代数：%.9g\t',meangen);
fprintf(file,'获得最优适应度的平均时间：%.9g\t',avg_besttime);
fprintf(file,'获得最优适应度的平均评估次数：%d\n\n', avg_evaluate_num);
fprintf(file,'平均测试时间：%gs\n\n',  meanruntime);
fprintf(file,'测试成功率：%g\n\n',  avg_success_count);
fprintf(file,'------------每一代平均最优适应度--------------------\n');
   for i=1:G
        fprintf(file,'%.9g\n',  avgBestPerGen(i));
    end
  fclose(file);


%根据轮盘赌算法选择新种群
%pop 种群
%fit 种群的适应度
%NP 种群个数
%D 个体维度
function pop_new=roulette(pop,fit,NP,D)
    pop_new=zeros(NP,D);
    maxf=max(fit);
    roulette_p=power(fit-maxf,2);
%     roulette_p=1./fit;
    sum_f=sum(roulette_p);
    roulette_p=roulette_p./sum_f;
    roulette_p=cumsum(roulette_p);
    for i=1:NP
        p=rand;
        for j=1:NP
            if p<roulette_p(j)
                pop_new(i,:)=pop(j,:);
                break;
            end
        end
    end

end





%变异
function pop_mutate=mutate(pop,NP,Pm,city_num)
    for i=1:NP
        r=rand;
        if r<Pm
             mutate_num=round(city_num*Pm);
             mutate_pos=randperm(city_num);
             for j=1:1  %每次只变异一个，若想变异多个则用mutate_num
                 temp=pop(i,mutate_pos(2*j));
                 pop(i,mutate_pos(2*j))=pop(i,mutate_pos(2*j-1));
                 pop(i,mutate_pos(2*j-1))=temp;
             end
        end
    end
    pop_mutate=pop;
end




%种群交叉，杂交方式采用了君主杂交，每个个体选连续几个分量与君主个体进行杂交，但同时也要保持杂交后路径序列不能出现重复
%pop 种群
%NP种群规模
%Pc杂交概率
%city_num 个体维度,这里是城市个数
function pop_cross=cross(pop,NP,Pc,city_num,bestPop)
    pop_cross=pop;
    for i=1:2:NP
%         pop_cross(i,:)=bestPop;%注释之后就变成单点交叉了
        cross_num=floor(Pc*city_num);%杂交的维数
        pos=unidrnd(city_num-cross_num+1);%从pos位置开始杂交，到pos+cross_num-1 
        common_individual_cross_part=pop_cross(i+1,pos:pos+cross_num-1);%普通个体用来杂交的部分
        emper_individual_cross_part=pop_cross(i,pos:pos+cross_num-1);%君主个体用来杂交的部分
        ind_cross=get_cross_ind(pop_cross(i+1,:),emper_individual_cross_part,cross_num,pos);%获取杂交后的普通个体
        emper_cross=get_cross_ind(pop_cross(i,:),common_individual_cross_part,cross_num,pos);%获取杂交后的君主个体
        pop_cross(i+1,:)=ind_cross;
        pop_cross(i,:)=emper_cross;
    end
   
end

%返回杂交后的个体
%individual 用于杂交的个体
%cross_part 用于和杂交个体进行杂交的部分，即替换杂交个体中的对应位置的分量
%cross_num 杂交分量的个数
%pos 杂交的起始位置
function ind_cross=get_cross_ind(individual,cross_part,cross_num,pos)
    city_num=length(individual);       
    special_flag=zeros(1,city_num);
        %标记与君主杂交的个体中需要被替换的分量，即替换后会缺失的分量
        special_flag(1,pos:pos+cross_num-1)=1;
        for k=1:cross_num
            index=find(individual==cross_part(k));
            if special_flag(index(1,1))==1%如果替换后会发生重复的分量的位置在待替换的位置上，则说明替换后该分量不会重复，同时因为也不会缺失因此置为0
                special_flag(index(1,1))=0;%置为0表示，替换后，该分量即不重复也不会缺失
            else
                special_flag(index(1,1))=2;%用2标记的分量，表示杂交后，个体中用2标记的分量表示会和杂交替换过来的分量重复
            end
        end
        %special_flag数组中，值为1与值为2的元素个数是一致的
        index1=find(special_flag==1);%1表示杂交后会缺失的分量
        index2=find(special_flag==2);%2表示杂交后会重复的分量
        len=length(index1);
        for m=1:len
            %用标记1的分量替换掉标记2的分量这样杂交后的个体就不会出现重复的分量
            individual(index2(1,m))=individual(index1(1,m));
        end
        individual(1,pos:pos+cross_num-1)=cross_part;
        ind_cross=individual;
end


%种群初始化
%NP 种群规模
%city_num 城市数量
function pop=initialPop(NP,city_num)
    pop=zeros(NP,city_num);
    for i=1:NP
        pop(i,:)=randperm(city_num);
    end
end
%适应度评价
%pop 记录了每个个体的路径
%NP 种群规模
%Dist 记录任意两个城市的距离
function fit=fitnessEvaluation(pop,NP,Dist)
    fit=zeros(NP,1);
    for i=1:NP
        fit(i)=getFit(pop(i,:),Dist);
    end
end

%返回一个个体的路径长度
%Path 记录一个个体的遍历路径
%D 所有城市任意两个顶点的距离
function fit=getFit(Path,D)
    len=length(Path);
    fit=0;
    for i=2:len
        dist=D(Path(i-1),Path(i));
        fit=fit+dist;
    end
    fit=fit+D(Path(1),Path(len));
end


%返回任意两个城市的距离矩阵
%n 城市的个数
%C 记录了各个城市的坐标
function Dist=getDistance(C,n)
    D=zeros(n,n);%用来记录任意两个城市的距离
    for i=1:n
        for j=i:n
            if i==j
                D(i,j)=0;
            else
                D(i,j)=sqrt((C(i,1)-C(j,1))^2+(C(i,2)-C(j,2))^2);
                D(j,i)=D(i,j);
            end
        end
    end
    Dist=D;
end


