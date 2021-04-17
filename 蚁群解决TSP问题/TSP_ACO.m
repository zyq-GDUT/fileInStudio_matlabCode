close all;
clc;
%信息素：信息素与启发式因子用于求蚂蚁选择下一个城市的概率
%启发式因子：这里用的是两个城市之间距离的倒数
%蚁群算法一次迭代的流程：m个蚂蚁随机落在n个城市中。每个蚂蚁的目的就是遍历完所有城市，且只经过各个城市一次。
%因此蚂蚁初始化随机落在一个城市，就是这个蚂蚁的起始遍历点。
%每个蚂蚁根据当前所在的城市与还没遍历的城市都能通过信息素与启发式因子求得一个概率，
%也就是一个蚂蚁目前所在一个城市，然后都是根据概率来选择下一个城市。选不同城市概率是不一样的。
%各个蚂蚁在寻找路径时是互不影响的。等到所有蚂蚁都遍历完城市，一次迭代算是完成。
%然后保存代的最优路径，同时计算每只蚂蚁所遍历路径的信息素增量。也就是，一个蚂蚁遍历的路径中，是分为一段一段的。
%然后根据这只蚂蚁遍历路径的长度给每段路径添加一个特殊值，这个值就是该段路的信息素增量。可以用一个二维矩阵来保存任意一段路的信息素增量。
%直到当前迭代的每个蚂蚁都计算完信息素增量，则更新信息素的二维矩阵，所有信息素先减少一定值（该值由信息素蒸发系数决定）再加上信息素增量，即是更新完信息素矩阵

%蚁群算法主要是解决组合优化问题
%解空间任意两个元素的关系用信息素和启发式因子来描述
%解空间的两个元素的信息素与启发式因子用于生成这两个元素之间的选择概率，即若当前已经将元素i作为解了，因为元素i与其他元素都会存在一个选择概率，
%概率越大，就越优可能成为i之后的一个解元素。
%任意两个城市的启发式因子在迭代过程中是不会变的。
%任意两个城市的信息素会随着迭代次数的增加而改变，从而改变每两个城市的选择概率
m=50;%蚂蚁个数
alpha=1;%信息素的权重
beta=5;%启发式因子的权重，这里是距离倒数的权重
Q=100;%用于求信息素增量
rho=0.1;%信息素蒸发系数，即每次迭代后当前代的信息素会降低百分之10
evaluateCount=10^5;%最大迭代次数
G=ceil(evaluateCount/m);%迭代次数
wonderful=15609.5;
error=10^-6;
MAXRUNS=20;
run=0;
%各个城市的坐标，用于计算相互之间的距离%
C=[1304 2312;3639 1315;4177 2244;3712 1399;3488 1535;3326 1556;...
3238 1229;4196 1044;4312  790;4386  570;3007 1970;2562 1756;...
2788 1491;2381 1676;1332  695;3715 1678;3918 2179;4061 2370;...
3780 2212;3676 2578;4029 2838;4263 2931;3429 1908;3507 2376;...
3394 2643;3439 3201;2935 3240;3140 3550;2545 2357;2778 2826;...
2370 2975];                
n=size(C,1);
D=getDistance(C,n);%获取任意两个城市的距离%
avgBestPerGen=zeros(G,1);%用来记录每次测试的平均每代最优
singleTestTime=zeros(MAXRUNS,1);%用来记录每次测试的时间
getBestTime=zeros(MAXRUNS,1);%用来记录每次测试中获取到最优值的时间
BestFitPerTest=zeros(MAXRUNS,1);%用来记录每次测试的最优适应度
genNum=zeros(MAXRUNS,1);%用来记录每次测试的最优代数
while run<MAXRUNS 
    Eta=1./D;%任意两个城市的启发式因子为距离的倒数
    Tau=ones(n,n);%信息素矩阵为n*n的矩阵，每个元素初始值为1，即每段路的信息素初始值为1，（i,j）元素值表示第i个元素到第j个元素的初始信息素
    R_best=zeros(G,n);%用来保存每代的最优路径
    Fit_best=zeros(G,1);
    ant_pop=zeros(m,n);%用来保存每次迭代蚂蚁遍历的路径
    ec=0;
    gen=1;
    % while gen<=G
    bestTimePerGen=zeros(G,1);%每代获取到最优的时间
    time1=clock;
    while ec<evaluateCount
        %让m只蚂蚁随机落在n个城市中，即初始化每只蚂蚁的第一个城市%
        ant_pop(:,1)=getrandpos(m,n);
        for i=2:n
            for j=1:m
                %找出还没被访问的城市
                visited=ant_pop(j,1:i-1);
                no_visited_city=get_novisited_city(visited,i-1,n);
                %计算当前城市到其他城市的概率
                PToOtherCity=getPToOtherCity(no_visited_city,ant_pop(j,i-1),Tau,Eta,alpha,beta);
                %构造轮盘赌选择器
                Pcum=cumsum(PToOtherCity);
                select=find(rand<=Pcum);
                ant_pop(j,i)=no_visited_city(select(1));  
            end
        end
        %这次迭代的保存最优值，与最优解
        [bestfit,bestPath,antPopFit]=getBest(ant_pop,D);
        ec=ec+m;
        R_best(gen,:)=bestPath;
        Fit_best(gen,1)=bestfit;
        if gen>1
            if Fit_best(gen,1)>Fit_best(gen-1,1)
                Fit_best(gen,1)=Fit_best(gen-1,1);
                R_best(gen,:)=R_best(gen-1,:);
                ant_pop(1,:)=R_best(gen,:);
            end
        end
        getbestTime=clock;
        bestTimePerGen(gen,1)=etime(getbestTime,time1);
        fprintf("第%d代最短路径为：%g\n",gen,Fit_best(gen,1));
%         drawFigure(n,R_best(gen,:),C,Fit_best(gen,1));
        gen=gen+1;
        %计算信息素增量
        delta_tau=getDeltaTau(ant_pop,antPopFit,Q,n);
        %更新信息素
        Tau=update_tau(delta_tau,Tau,rho);
    end
    time2=clock;
    [best,index]=min(Fit_best);
    getBestTime(run+1,1)=bestTimePerGen(index,1);
    BestFitPerTest(run+1,1)=best;
    genNum(run+1,1)=index;
    singleTestTime(run+1,1)=etime(time2,time1);
    avgBestPerGen=avgBestPerGen+Fit_best;
    run=run+1;
end
avgBestPerGen= avgBestPerGen/MAXRUNS;
meanfit=mean(BestFitPerTest);
std_bestfit=std(BestFitPerTest);
bestfit=min(BestFitPerTest);%所有测试中最优值的最优值
worstfit=max(BestFitPerTest);%所有测试中最优值中的最差值
meangen=mean(genNum);%每次测试首次获得最优适应度的代数
meanruntime=mean(singleTestTime);%单次测试的平均时间
avg_evaluate_num=meangen*m;
avg_besttime=mean(getBestTime);
success_count=0;
for i=1:MAXRUNS
    if (BestFitPerTest(i)-wonderful)<=error
        success_count=success_count+1;
    end
end
avg_success_count=success_count/MAXRUNS;



file=fopen('testdata.txt','w');
fprintf(file,'算法名：ACO算法解决TSP问题\n\n测试次数：%d\t迭代次数：%d\t总的评价次数：%d\n\n蚂蚁个数：%d\t问题维度：%d\n\n信息素权重alpha：%.3f\t启发式因子权重beta:%.3f\n\n蒸发系数rho:%.3f\t信息素增量权重Q：%.3f\n\n',MAXRUNS,G,evaluateCount,m,n,alpha,beta,rho,Q);
fprintf(file,'平均最优适应度：%.9g\t',meanfit);
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
        fprintf(file,'%.9g\n',avgBestPerGen(i));
    end
  fclose(file);
 
 
 
 
 
%画出一条路径
%len是路径的长度
%path 记录了整条路径的编号
%C 记录了每个点的坐标
%fitbest 该路径的距离
function drawFigure(len,path,C,fitbest)
    for i=1:len-1
        plot([C(path(i),1),C(path(i+1),1)],[C(path(i),2),C(path(i+1),2)],'bo-');
        hold on;%让多段线共存，目的是绘画出整条路径
    end
    plot([C(path(len),1),C(path(1),1)],[C(path(len),2),C(path(1),2)],'ro-');
    hold off;%完整的路径绘画完毕后，就可以让下一代的图覆盖
    title(['该路径的距离是：',num2str(fitbest)]);
    pause(0.005);
end


%更新信息素矩阵
%tau 信息素矩阵
%delta_tau 信息素增量矩阵
%rho 信息素蒸发系数 
function uptodate_tau=update_tau(delta_tau,tau,rho)
    uptodate_tau=tau*(1-rho)+delta_tau;
end


%计算这次迭代后，任意两个城市的信息素增量
%ant_pop 记录了每只蚂蚁的遍历路径
%antPopFit 记录了每只蚂蚁遍历路径的长度
%Q 信息素增量的权重
%n 城市的数量
function delta_tau=getDeltaTau(ant_pop,antPopFit,Q,n)
    delta_tau=zeros(n,n);
    m=size(ant_pop,1);
    for i=1:m
        for j=2:n
            index1=ant_pop(i,j-1);
            index2=ant_pop(i,j);
            delta_tau(index1,index2)=delta_tau(index1,index2)+Q/antPopFit(i);
            %从1到2信息素浓度高，说明1到2大概率是距离短的，因此2到1距离也是短的，因此返过来2到1的信息素浓度也应该高
%             delta_tau(index2,index1)=delta_tau(index2,index1)+Q/antPopFit(i);
        end
        delta_tau(ant_pop(i,1),ant_pop(i,n))=delta_tau(ant_pop(i,1),ant_pop(i,n))+Q/antPopFit(i);
%         delta_tau(ant_pop(i,n),ant_pop(i,1))=delta_tau(ant_pop(i,n),ant_pop(i,1))+Q/antPopFit(i);
    end
    
end


%返回当代最小的路径长度，对应的路径,还有各个蚂蚁的路径长度
function [bestfit,bestPath,antPopFit]=getBest(ant_pop,D)
   
    antPopFit=getAntPopFit(ant_pop,D);
    num=size(antPopFit,1);
    min_fit=antPopFit(1,1);
    min_index=1;
    for i=1:num
        if(antPopFit(i,1)<min_fit)
            min_fit=antPopFit(i,1);
            min_index=i;
        end
    end
    bestPath=ant_pop(min_index,:);
    bestfit=min_fit;
end

%返回所有蚂蚁的路径长度
function antPopFit=getAntPopFit(ant_pop,D)
    m=size(ant_pop,1);
    antPopFit=zeros(m,1);
    for i=1:m
        antPopFit(i,1)=getFit(ant_pop(i,:),D);
    end
end



%返回一个蚂蚁遍历路径的路径长度
%antPath 记录一只蚂蚁的遍历路径
%D 所有城市任意两个顶点的距离
function fit=getFit(antPath,D)
    len=length(antPath);
    fit=0;
    for i=2:len
        dist=D(antPath(i-1),antPath(i));
        fit=fit+dist;
    end
    fit=fit+D(antPath(1),antPath(len));
end


%返回当前城市到其他城市的概率矩阵
%noVisitedCity 记录了还没被访问的城市
%currentCity 当前所在的城市
%Eta 记录了任意两个城市的启发式因子
%Tau 记录了当前任意两个城市的信息素
%alpha 信息素的权重
%beta 启发式因子的权重
function PToOtherCity=getPToOtherCity(noVisitedCity,currentCity,Tau,Eta,alpha,beta)
    len=length(noVisitedCity);
    PToOtherCity=zeros(1,len);
    for i=1:len
        PToOtherCity(1,i)=(Tau(currentCity,noVisitedCity(1,i))^alpha)*(Eta(currentCity,noVisitedCity(1,i))^beta);
    end
  PToOtherCity=PToOtherCity/sum(PToOtherCity);
end


%用于返回待访问的城市列表
%visitedcity 已访问的城市
%n visitedcity的元素个数，即已访问城市的个数
%allcitynum 1~allcitynum为每个城市的编号

function citylist=get_novisited_city(visitedCity,visitedNum,allCityNum)
    j=1;        
    citylist=zeros(1,allCityNum-visitedNum);
    for i=1:allCityNum
        len=length(find(visitedCity==i));
        if len==0
            citylist(j)=i;
            j=j+1;
        end
    end
    
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






 
%返回每个蚂蚁的初始城市编号
%ant_num 蚂蚁的数量
%city_num 城市的数量
function result=getrandpos(ant_num,city_num)
 n=ceil(ant_num/city_num);%向上取整，即求成蚂蚁数量是城市数量的n倍
 randpos=zeros(1,n*city_num);
 for i=1:n
%     randpos=[randpos,randperm(city_num)];
        randpos(1,(i-1)*city_num+1:i*city_num)=randperm(city_num);
 end
 result=randpos(1,1:ant_num)';
end





