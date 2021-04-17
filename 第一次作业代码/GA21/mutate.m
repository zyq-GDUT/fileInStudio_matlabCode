%±äÒì²Ù×÷
%
function pop_mutate=mutate(pop,NP,Pm,D,lb,ub)
    mutate_num=round(NP*Pm);
    for i=1:mutate_num
        rp=randi([1,NP],1,1);
        for j=1:(D*Pm)
            rd=randi([1,D],1,1);
            pop(rp,rd)=lb+(ub-lb)*rand;
        end
    end
    pop_mutate=pop;
end