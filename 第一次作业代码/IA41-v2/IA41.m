%%%%%初始化%%%%%%
clear all;
clear all;
clc;
D=10; %解的维度，也就是抗体的维度
NP=100; %抗体的个数
ub=20;  %定义域上限
lb=-20; %定义域下限
G=500; %最大免疫代数
pm=0.7;  %变异的概率

global alfa;
alfa=1;%激励度的亲和度系数
global belta;
belta=-1; %激励度的浓度系数
gen=0;  %代数
Nc1=10; %克隆个数
global detas;
detas=0.1;
global deta0;
deta0=1*ub; %领域范围初值
bestfes=10^5;%最大评价次数 
fes=0;%当前评价次数

FIT=zeros(NP,1); %返回一个NP*1的全零矩阵
%%初始化种群%%%
f=rand(D,NP)*(ub-lb)+lb;
for np=1:NP
 FIT(np)=func1(f(:,np)); %计算抗体的亲和度，即适应度
end
fes=fes+NP;
%计算抗体浓度%
ND=get_consistence(f,NP);
 bestfit_perGen=zeros(G,1);
%用抗体的亲和度和浓度计算激励度%
FIT=alfa*FIT-belta*ND;
%激励度按升序排列%
[SortFIT,Index]=sort(FIT);
Sortf=f(:,Index);%得到的sortf种群与激励度数组SortFIT一一对应
aFIT=zeros(NP/2,1);
af=zeros(D,NP/2);
bFIT=zeros(NP/2,1);
%while gen<G
while fes<bestfes
    for i=1:NP/2 %变异选择即取前一半的个体，即激励度最优的一半个体。
        a=Sortf(:,i);
        Na=repmat(a,1,Nc1);%按列复制Nc1个抗体a
        deta=deta0/(gen+1); %用于局部搜索的领域会随着迭代次数的增加而变小
        %每个克隆的抗体的每一维度都按概率进行变异%
        for ii=1:Nc1
            for j=1:D
                if rand<pm
                    Na(j,ii)= Na(j,ii)+(rand-0.5)*deta;
                end
            end
        end
         Na(:,1)=Sortf(:,i);%保证了局部搜索后的抗体>=原克隆抗体
        %%%克隆抑制，即找出10个克隆抗体中适应度最优的抗体%%%
        NaFIT=zeros(Nc1,1);
        for j=1:Nc1
            NaFIT(j)=func1(Na(:,j));
        end  
        fes=fes+Nc1;
            [NaSortFIT,Index]=sort(NaFIT);
            aFIT(i)=NaSortFIT(1);
            NaSortf=Na(:,Index);
            af(:,i)=NaSortf(:,1);
    end
    %计算免疫种群的激励度%
    %计算免疫种群的浓度
    afND=get_consistence(af,NP/2);    
    afJLD=alfa*aFIT-belta*afND;%免疫个体的激励度%
    %af是免疫抗体的信息，aFIT是免疫抗体对应的亲和度，afJLD是免疫抗体对应的激励度%
    %种群刷新，重新随机生成剩余的一半抗体%
    bf=rand(D,NP/2)*(ub-lb)+lb;
    for np=1:NP/2
        bFIT(np)=func1(bf(:,np));
    end
    fes=fes+(NP/2);
    bfND=get_consistence(bf,NP/2); %后一半抗体bf的浓度
    bfJLD=alfa*bFIT-belta*bfND;
    %免疫种群与新生种群合并%
    f1=[af,bf];%按列合并
    FIT1=[aFIT;bFIT];%按行合并
    JLD=[afJLD;bfJLD];%按行合并
    [bestfit,bestindex]=min(FIT1,[],1);%按行找最小值%
    gen=gen+1;
    bestfit_perGen(gen)=bestfit;
    best_individual=f1(:,bestindex);
    [SortJLD,Index]=sort(JLD);
    Sortf=f1(:,Index);
    
end
 for i=1:G
 fprintf("第%d代最优：%g\n", i,bestfit_perGen(i));
 end






        
