%适应度评估
%fit保存了每个个体的适应度
function fit=evaluate(pop,NP)
    fit = zeros(1,NP);
    for i=1:NP    
        fit(i)=evaluationfun(pop(i,:));
    end
end