%根据轮盘赌算法选择新种群
%pop 种群
%fit 种群的适应度
%NP 种群个数
%D 个体维度
function pop_new=roulette(pop,fit,NP,D)
    pop_new=zeros(NP,D);
    maxf=max(fit);
    roulette_p=power(fit-maxf,2);
%     roulette_p=1./fit;
    sum_f=sum(roulette_p);
    roulette_p=roulette_p./sum_f;
    roulette_p=cumsum(roulette_p);
    for i=1:NP
        p=rand;
        for j=1:NP
            if p<roulette_p(j)
                pop_new(i,:)=pop(j,:);
                break;
            end
        end
    end

end