%求每个个体的浓度
%pop是种群
%N是种群规模
function ND=get_consistence(pop,NP)
    global detas;
    ND=zeros(NP,1);
    nd=zeros(NP,1);
    for np=1:NP
        for j=1:NP
           nd(j)=sqrt(sum((pop(:,np)-pop(:,j)).^2));
           if nd(j)<detas
               nd(j)=1;
           else
               nd(j)=0;
           end
        end
        ND(np)=sum(nd)/NP;
    end
end
        
