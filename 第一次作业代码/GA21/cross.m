%种群交叉
%pop 种群
%NP种群规模
%Pc杂交概率
%D 个体维度
function pop_cross=cross(pop,NP,Pc,D)
    for i=1:2:NP
        p=rand;
        if p<Pc
            p2=randi([0,1],1,D);
            for j=1:D
                if p2(j)==1
                 temp=pop(i+1,j);
                 pop(i+1,j)=pop(i,j);
                 pop(i,j)=temp;
                end
            end
        end
    end
    pop_cross=pop;

end