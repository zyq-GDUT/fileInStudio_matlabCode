
%% 
clear all;
close all;
clc;
global NP D Xx Xs G

D=10;
NP=100;
Xs=20;
Xx=-20;
G=2000;
f=zeros(D,NP);
nf=zeros(D,NP);
Pc=0.8;
Pm=0.1;
%% ������ݵĶ��� %%
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
%% ��ʼ���� %%
for run =1:MAXRUN
    fes = 0;
    t0 = clock;
    %% ��ʼ����Ⱥ %%
    f=rand(D,NP)*(Xs-Xx)+Xx; %��Ⱥ��ʼ��
    %% ������Ӧ�� %%
    FIT = zeros(1,NP);
    for np=1:NP
        FIT(np) = funcl(f(:,np));
        fes=fes+1;
    end
    [SortFIT,Index]=sort(FIT);
    Sortf = f(:,Index);
    %% �������׼�� %%
    gbestFIT = SortFIT(1);  %���ѡ��һ�������Ŀ�꺯��ֵ�䵱��ʼ���������Ž�ֵ
    bestFEs = fes;
    bestGen = 0;
    besttime = 0;
    onerunfile = fopen(['F1_run' num2str(run) '.txt'],'w'); %w ��Ϊa��ʾ���   
    tall = 0;        %���β�������ʱ
    gen=0;           %���ߴ���
    flagOK = false;  %��������ڵ�����ֵʱתΪtrue
%%%%%%%%%%%%%%%%%%%�Ŵ��㷨ѭ��%%%%%%%%%%%%%%%%%%%%
     tt1=clock;
    for gen=1:G     
         nf = crossover(Sortf,Pc);  
         nf = mutation(nf,Pm);
         [Sortf,SortFIT]=selection(Sortf,nf);
         trace(gen) = SortFIT(1);
         if trace(gen) < gbestFIT             %�Ƚ��Ƿ��ҵ����ŵ��������Ž�+++++
                gbestFIT = trace(gen);           %�����������Ž�ֵ
                bestGen = gen;                   %�״λ���������Ž�ֵ��Ӧ��ƽ������
                bestFEs = fes;                     %�״λ���������Ž�ֵ��Ӧ�ĺ������۴���
                besttime = tall + etime(clock,t0); %��Ӧ����ʱ
                if(~flagOK && gbestFIT < error)    %�״��ҵ����㾫��error�ڵĽ�
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
     timecount(run)=etime(tt2,tt1);
     fclose(onerunfile);
        tallrecord = [tallrecord,tall];
        bestfrecord = [bestfrecord, trace(gen)];
        funfile = fopen('F1total.txt','a'); %w ��Ϊa��ʾ���
        fprintf('%d\t%g\t%g\r\n',run,trace(gen),tall)
        fprintf(funfile,'%d\t%g\t%g\t%g\t%g\t%g\r\n',run,trace(gen),bestGen,bestFEs,besttime,tall);
        fclose(funfile);
        meantrace = [meantrace; trace];  %���кϲ����������м�¼
        bestGenRecord = [bestGenRecord,bestGen];
        bestFEsRecord = [bestFEsRecord,bestFEs];
        bestTimeRecord = [bestTimeRecord,besttime];
end
timecount
meantrace = mean(meantrace);     %����ȡ��ֵ�����ƽ����������
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
        
%Bestf = Sortf(:,1);
%trace(end);
%figure
%plot(trace)
%xlabel('��������')
%ylabel('Ŀ�꺯��ֵ')
%title('��Ӧ�Ƚ�������')
     
  