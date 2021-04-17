function cpop = crossover(pop, Pc)
    global D NP
     Emper=pop(:,1);  
     NoPoint = round(D*Pc);  %交叉个数
     PoPoint = randi([1 D],NoPoint,NP/2); 
     nf= pop;
     %君主交叉
     for i = 1:NP/2
         nf(:,2*i-1)= Emper;
         nf(:,2*i) = pop(:,2*i);
         for k =1:NoPoint
            nf(PoPoint(k,i),2*i-1)=nf(PoPoint(k,i),2*i);
            nf(PoPoint(k,i),2*i)=Emper(PoPoint(k,i));
         end
     end
     cpop = nf;
end

