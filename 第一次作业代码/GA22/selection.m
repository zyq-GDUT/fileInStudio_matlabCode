function [npop,nfit] = selection(pop, mpop)

     global NP fes;
     for np = 1:NP
         NFIT(np)=funcl(mpop(:,np));
     end
     for np = 1:NP
         SortFIT(np)=funcl(pop(:,np));
         fes=fes+1;
     end
     [NSortFIT,Index] = sort(NFIT);
     NSortf = mpop(:,Index);
     %产生新种群
     f1 = [pop,NSortf];
     FIT1 = [SortFIT,NSortFIT];
     [SortFIT1,Index] = sort(FIT1);
     Sortf1 = f1(:,Index);
     SortFIT = SortFIT1(1:NP);
     pop = Sortf1(:,1:NP);
     npop = pop;
     nfit = SortFIT;
end

