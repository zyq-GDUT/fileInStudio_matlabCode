clear all;
close all;
clc;
%%初始化参数
NP=100;
D=10;
Pc=0.8;
Pm=0.1;
G=2000;%迭代次数
ub=20;%上界
lb=-20;%下界

%% 输出数据的定义 %%
error=10^-6;
MAXRUN =1;
bestFitEveryGen=zeros(MAXRUN,G);   
global evaluateNum;
evaluateNum=zeros(MAXRUN);
time=zeros(1,MAXRUN);
besttime=zeros(MAXRUN,G);
% CoreProcess=coreProcess;
for run=1:MAXRUN  
    t1=clock;
    pop=initpop(NP,D,lb,ub);
    for i=1:G
        fit=evaluate(pop,NP);
%         evaluateNum(run)=evaluateNum(run)+NP;%总的评价次数%
        fit_min=min(fit);
        bestt=clock;
        besttime(run,i)=etime(bestt,t1);
        index=find(fit_min==fit);
        best_individual=pop(index(1,1),:);
        fprintf("当前最优个体：");
        best_individual
        bestFitEveryGen(run,i)=fit_min;%保存每代最优适应度%
        pop=roulette(pop,fit,NP,D);
        pop=cross(pop,NP,Pc,D);
        pop=mutate(pop,NP,Pm,D,lb,ub);
        pop(1,:)=best_individual;
    end
    t2=clock;
    time(run)=etime(t2,t1);
    fprintf('第%d次测试：最优值为=%.3f，用时%.3f\n',run,bestFitEveryGen(run,G), time(run));

end
% time
  [minfit,index]=min(bestFitEveryGen,[],2); 
   meanfit=mean(minfit(:));%每次测试的最优适应度
   meangen=mean(index(:));%每次测试首次获得最优适应度的代数
   meanruntime=mean(time(:));%单次测试的平均时间
   bestfit=min(minfit);%所有测试中最优值的最优值
   worstfit=max(minfit);%所有测试中最优值中的最差值
   std_bestfit=std(minfit);%标准差
   sum_besttime=0;sum_evaluate=0;success_count=0;
   for i=1:MAXRUN
    sum_besttime=sum_besttime+besttime(i,index(i,1));
    sum_evaluate=sum_evaluate+index(i,1)*NP;
    if minfit(i,1)<=error
        success_count=success_count+1;
    end
   end
   avg_besttime=sum_besttime/MAXRUN;
   avg_evaluate_num=sum_evaluate/MAXRUN;
   avg_success_count=success_count/MAXRUN;
%    meanfit
%    std_bestfit
%    bestfit
%    worstfit
%    meangen
%    meanruntime
%    avg_besttime
%    avg_evaluate_num
%    avg_success_count
   
   file=fopen('testdata.txt','w');
   fprintf(file,'算法名：GA21\n\n测试次数：%d\t迭代次数：%d\n\n种群规模：%d\t交叉概率：%.3f\t变异概率:%.3f\n\n',MAXRUN,G,NP,Pc,Pm);
   fprintf(file,'平均最优适应度：%.3f\t',meanfit);
   fprintf(file,'最优适应度的标准差：%.3f\n\n',std_bestfit);
   fprintf(file,'最优适应度的最优值：%.3f\t',bestfit);
   fprintf(file,'最优适应度的最差值：%.3f\n\n',worstfit);
   fprintf(file,'获得最优适应度的平均代数：%.3f\t',meangen);
   fprintf(file,'获得最优适应度的平均时间：%.3f\t',avg_besttime);
   fprintf(file,'获得最优适应度的平均评估次数：%d\n\n', avg_evaluate_num);
   fprintf(file,'平均测试时间：%.3fs\n\n',  meanruntime);
   fprintf(file,'测试成功率：%.3f\n\n',  avg_success_count);
   fprintf(file,'------------每一代平均最优适应度--------------------\n');
   avgbestfit_pergGen=mean(bestFitEveryGen,1);%按列平均
 
    for i=1:G
        fprintf(file,'第%d代的平均最优适应度：%.3g\n',  i,avgbestfit_pergGen(i));
    end
   fclose(file);
