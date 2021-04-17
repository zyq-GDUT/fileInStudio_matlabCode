%% ��ʼ�� %%
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
    tt1=clock;
    fes = 0;
    t0 = clock;
    %% ��ʼ����Ⱥ %%
    x=zeros(D,NP);
    v=zeros(D,NP);
    u=zeros(D,NP);
    x=rand(D,NP)*(Xs-Xx)+Xx;
    %% ������Ӧ�� %%
    for m =1:NP
        Ob(m)=funcl(x(:,m));
        fes=fes+1;
    end
    trace(1)=min(Ob);
    %% �������׼�� %%
    gbestFIT = Ob(1);  %���ѡ��һ�������Ŀ�꺯��ֵ�䵱��ʼ���������Ž�ֵ
    bestGen = 0;
    bestFEs = fes;
    besttime = 0;
    onerunfile = fopen(['F1_run' num2str(run) '.txt'],'w'); %w ��Ϊa��ʾ���   
    tall = 0;        %���β�������ʱ
    gen=0;           %���ߴ���
    flagOK = false;  %��������ڵ�����ֵʱתΪtrue
    %% ���ѭ�� %%
    for gen=1:G
        v = mutation32(x,NP,F);
        u = crosser32(D,CR,v,x,NP,Xx,Xs);
        selection32(NP,u,gen)
        if trace(gen+1) < gbestFIT             %�Ƚ��Ƿ��ҵ����ŵ��������Ž�+++++
            gbestFIT = trace(gen+1);           %�����������Ž�ֵ
            bestGen = gen+1;                   %�״λ���������Ž�ֵ��Ӧ��ƽ������
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
    timecount(run)=etime(tt2,tt1);
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
%% ��ͼ %%
%figure
%plot(trace);
%xlabel('��������')
%ylabel('Ŀ�꺯��ֵ')
%title('DEĿ�꺯������')

