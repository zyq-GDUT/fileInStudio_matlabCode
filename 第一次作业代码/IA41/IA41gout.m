%%%%%%%%%%%%%%%%%�����㷨������ֵ%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%��ʼ��%%%%%%%%%%%%%%%%%%%%%%%%%
clear all;                                %������б���
close all;                                %��ͼ
%clc;                                      %����
global fes
MAXRUN = 20;                %�������Դ���
D=10;                                     %���߸���ά��
NP=100;                                   %���߸�����Ŀ
Xs=20;                                    %ȡֵ����
Xx=-20;                                   %ȡֵ����
G=1000;                                    %������ߴ���
pm=0.7;                                   %�������
alfa=1;                                   %������ϵ��
belta=1;                                  %������ϵ��
detas=0.2;                                %���ƶ���ֵ
Ncl=10;                                   %��¡����
deta0=1*Xs;                               %����Χ��ֵ
error = 10^-6;                  %��ý�����ֵ

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
    %%%%%%%%%%%%%%%%%%%%%%%��ʼ��Ⱥ%%%%%%%%%%%%%%%%%%%%%%%%
    f=rand(D,NP)*(Xs-Xx)+Xx;            %NP������������
    for np=1:NP
        FIT(np)=func1(f(:,np));         %�׺Ͷ� = Ŀ�꺯��ֵ
    end
    gbestFIT = FIT(1);      %���ѡ��һ�������Ŀ�꺯��ֵ�䵱��ʼ���������Ž�ֵ +++++++++
    gbextX = f(:,1);        %���ѡ��һ��������Ϊ�������Ž�+++++++++++
    bestGen = 0;
    bestFEs = fes;
    besttime = 0;
    %%%%%%%%%%%%%%%%%�������Ũ�Ⱥͼ�����%%%%%%%%%%%%%%%%%%%
    for np=1:NP
        for j=1:NP
            nd(j)=sqrt(sum((f(:,np)-f(:,j)).^2)); %������׺Ͷ�P63��4.3��
            if nd(j)<detas          %P62��4.2���е����ƶ���ֵ��=0.2��
                nd(j)=1;            %����i�͸���j����
            else
                nd(j)=0;            %����i�͸���j����Ƚ�Զ��������
            end
        end
        ND(np)=sum(nd)/NP;          %����Ũ��P62��4.1����Խ����Ũ��Խ��
    end
    FIT = alfa*FIT- belta*ND; %����Ŀ�꺯��ֵ�Ϳ���Ũ�ȵõ����弤����P63��4.6��
    %%%%%%%%%%%%%%%%%%%�����Ȱ���������%%%%%%%%%%%%%%%%%%%%%%
    [SortFIT,Index]=sort(FIT);
    Sortf=f(:,Index);               %����ֵԽСԽ��
    
    onerunfile = fopen(['F1_run' num2str(run) '.txt'],'w'); %w ��Ϊa��ʾ���   
    t0 = clock;
    tall = 0;                       %���β�������ʱ
    gen=0;                          %���ߴ���
    flagOK = false;
    %%%%%%%%%%%%%%%%%%%%%%%%����ѭ��%%%%%%%%%%%%%%%%%%%%%%%%
    while gen<G
        for i=1:NP/2    %ѡ������ǰNP/2����������߲���%%%%%%%%%%%
            a=Sortf(:,i);         %��i������ı���a
            Na=repmat(a,1,Ncl);   %��¡����Nc1�ݣ���1��Nc1�У�Nc1=10��
            deta=deta0/gen;       %��������Χֵ�����gen��С
            for j=1:Ncl           %�Կ�¡��ÿһ������
                for ii=1:D
                    %%%%%%%%%%%%%%%%%����%%%%%%%%%%%%%%%%%%%
                    if rand<pm    %�������������
                        Na(ii,j)=Na(ii,j)+(rand-0.5)*deta; %�����ڵľֲ�����P64��4.9��
                    end
                    %%%%%%%%%%%%%%�߽���������%%%%%%%%%%%%%%%
                    if (Na(ii,j)>Xs)  ||  (Na(ii,j)<Xx)
                        Na(ii,j)=rand * (Xs-Xx)+Xx;     %���ֵ�滻
                    end
                end
            end
            
            Na(:,1)=Sortf(:,i); %������¡Դ���壬�˴�������ԭ����ȡ����¡�ĵ�һ�����츱��
            
            %%%%%%%%%%��¡���ƣ������׺Ͷ���ߵĸ���%%%%%%%%%%
            for j=1:Ncl
                NaFIT(j)=func1(Na(:,j));    %��ÿ����¡��������һ����ԭ���壩��������
            end
            [NaSortFIT,Index]=sort(NaFIT);
            aFIT(i)=NaSortFIT(1);           %���ŵĺ�������ֵ������aFIT(i)
            NaSortf=Na(:,Index);
            af(:,i)=NaSortf(:,1);           %�����������Ÿ��屣����af�ĵ�i��
            
            if aFIT(i) < gbestFIT             %�Ƚ��Ƿ��ҵ����ŵ��������Ž�+++++
                gbestFIT = aFIT(i);           %�����������Ž�ֵ
                gbextX = af(:,i);          %�����������Ž�
                bestGen = gen;                 %�״λ���������Ž�ֵ��Ӧ��ƽ������
                bestFEs = fes;                 %�״λ���������Ž�ֵ��Ӧ�ĺ������۴���
                besttime = tall + etime(clock,t0);  %��Ӧ����ʱ
                if(~flagOK && gbestFIT < error) %�״��ҵ����㾫��error�ڵĽ�
                    flagOK = true;
                    okNum = okNum +1;
                end
            end
        end
        %%%%%%%%%%%%%%%%%%%%������Ⱥ������%%%%%%%%%%%%%%%%%%%
        for np=1:NP/2               %һ����Ⱥ
            for j=1:NP/2
                nda(j)=sqrt(sum((af(:,np)-af(:,j)).^2)); %������׺Ͷ�P63��4.3��
                if nda(j)<detas      %P62��4.2���е����ƶ���ֵ��=0.2��
                    nda(j)=1;
                else
                    nda(j)=0;
                end
            end
            aND(np)=sum(nda)/NP/2;  %�˴�Ϊһ����Ⱥ
        end
        aFIT =  alfa*aFIT-  belta*aND;  %�õ�һ����Ⱥ�Ŀ��弤����
        %%%%%%%%%%%%%%%%%%%%%%%��Ⱥˢ��%%%%%%%%%%%%%%%%%%%%%%%
        bf=rand(D,NP/2)*(Xs-Xx)+Xx;     %��һ����Ⱥ�������
        for np=1:NP/2
            bFIT(np)=func1(bf(:,np));
        end
        
        [ibestFIT,irank] = min(bFIT);      %���ε����к���ֵ���ŵ�ֵ+++++++++
        if ibestFIT < gbestFIT             %�Ƚ��Ƿ��ҵ����ŵ��������Ž�+++++
            gbestFIT = ibestFIT;           %�����������Ž�ֵ
            gbextX = bf(:,irank);          %�����������Ž�
            bestGen = gen;                 %�״λ���������Ž�ֵ��Ӧ��ƽ������
            bestFEs = fes;                 %�״λ���������Ž�ֵ��Ӧ�ĺ������۴���
            besttime = tall + etime(clock,t0);  %��Ӧ����ʱ
            if(~flagOK && gbestFIT < error) %�״��ҵ����㾫��error�ڵĽ�
                flagOK = true;
                okNum = okNum +1;
            end
        end
        
        %%%%%%%%%%%%%%%%%%%��������Ⱥ������%%%%%%%%%%%%%%%%%%%%
        for np=1:NP/2               %һ����Ⱥ
            for j=1:NP/2
                ndc(j)=sqrt(sum((bf(:,np)-bf(:,j)).^2));
                if ndc(j)<detas
                    ndc(j)=1;
                else
                    ndc(j)=0;
                end
            end
            bND(np)=sum(ndc)/NP/2;    %�˴�Ϊһ����Ⱥ
        end
        bFIT =  alfa*bFIT-  belta*bND; %�õ���һ����Ⱥ�Ŀ��弤����
        %%%%%%%%%%%%%%������Ⱥ��������Ⱥ�ϲ�%%%%%%%%%%%%%%%%%%%
        f1=[af,bf];                   %������Ⱥ�����кϲ�
        FIT1=[aFIT,bFIT];             %������Ⱥ���弤���Ⱥϲ�
        [SortFIT,Index]=sort(FIT1);   %������Ⱥ�����弤��������
        Sortf=f1(:,Index);            %�������������Ⱥ
        gen=gen+1;
        trace(gen)= gbestFIT; %func1(Sortf(:,1)); %��¼���ż����ȿ����Ӧ�ĺ�������ֵ++++
        
        tall = tall + etime(clock,t0);
        sprintf('%d\t%d\t%g\r\n',gen,fes,trace(gen)); %%%%%%out%%%%%%%%%%%
        fprintf(onerunfile,'%d\t%d\t%g\r\n',gen,fes,trace(gen));
        t0 = clock;
    end
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

meantrace = mean(meantrace);     %����ȡ��ֵ�����ƽ����������
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

%%%%%%%%%%%%%%%%%%%%%%%����Ż����%%%%%%%%%%%%%%%%%%%%%%%%
% Bestf=Sortf(:,1);          %���Ÿ������
% trace(end);                %����ֵ��end����ȡtrace�������һ����
% figure,plot(trace)
% xlabel('��������')
% ylabel('Ŀ�꺯��ֵ')
% title('�׺ͶȽ�������')