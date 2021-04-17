clear all;
close all;
clc;
%%��ʼ������
NP=100;
D=10;
Pc=0.8;
Pm=0.1;
G=2000;%��������
ub=20;%�Ͻ�
lb=-20;%�½�

%% ������ݵĶ��� %%
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
%         evaluateNum(run)=evaluateNum(run)+NP;%�ܵ����۴���%
        fit_min=min(fit);
        bestt=clock;
        besttime(run,i)=etime(bestt,t1);
        index=find(fit_min==fit);
        best_individual=pop(index(1,1),:);
        fprintf("��ǰ���Ÿ��壺");
        best_individual
        bestFitEveryGen(run,i)=fit_min;%����ÿ��������Ӧ��%
        pop=roulette(pop,fit,NP,D);
        pop=cross(pop,NP,Pc,D);
        pop=mutate(pop,NP,Pm,D,lb,ub);
        pop(1,:)=best_individual;
    end
    t2=clock;
    time(run)=etime(t2,t1);
    fprintf('��%d�β��ԣ�����ֵΪ=%.3f����ʱ%.3f\n',run,bestFitEveryGen(run,G), time(run));

end
% time
  [minfit,index]=min(bestFitEveryGen,[],2); 
   meanfit=mean(minfit(:));%ÿ�β��Ե�������Ӧ��
   meangen=mean(index(:));%ÿ�β����״λ��������Ӧ�ȵĴ���
   meanruntime=mean(time(:));%���β��Ե�ƽ��ʱ��
   bestfit=min(minfit);%���в���������ֵ������ֵ
   worstfit=max(minfit);%���в���������ֵ�е����ֵ
   std_bestfit=std(minfit);%��׼��
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
   fprintf(file,'�㷨����GA21\n\n���Դ�����%d\t����������%d\n\n��Ⱥ��ģ��%d\t������ʣ�%.3f\t�������:%.3f\n\n',MAXRUN,G,NP,Pc,Pm);
   fprintf(file,'ƽ��������Ӧ�ȣ�%.3f\t',meanfit);
   fprintf(file,'������Ӧ�ȵı�׼�%.3f\n\n',std_bestfit);
   fprintf(file,'������Ӧ�ȵ�����ֵ��%.3f\t',bestfit);
   fprintf(file,'������Ӧ�ȵ����ֵ��%.3f\n\n',worstfit);
   fprintf(file,'���������Ӧ�ȵ�ƽ��������%.3f\t',meangen);
   fprintf(file,'���������Ӧ�ȵ�ƽ��ʱ�䣺%.3f\t',avg_besttime);
   fprintf(file,'���������Ӧ�ȵ�ƽ������������%d\n\n', avg_evaluate_num);
   fprintf(file,'ƽ������ʱ�䣺%.3fs\n\n',  meanruntime);
   fprintf(file,'���Գɹ��ʣ�%.3f\n\n',  avg_success_count);
   fprintf(file,'------------ÿһ��ƽ��������Ӧ��--------------------\n');
   avgbestfit_pergGen=mean(bestFitEveryGen,1);%����ƽ��
 
    for i=1:G
        fprintf(file,'��%d����ƽ��������Ӧ�ȣ�%.3g\n',  i,avgbestfit_pergGen(i));
    end
   fclose(file);
