function u = crosser31(D,CR,v,x,NP,Xx,Xs)
%% 交叉 %%
    r= randi([1,D],1,1);
    for n=1:D
        cr=rand(1);
        if(cr<=CR)||(n==r)
            u(n,:)=v(n,:);
        else
            u(n,:)=x(n,:);
        end
    end
    %% 边界处理 %%
    for n =1:D
        for m=1:NP
            if(u(n,m)<Xx ||(u(n,m)>Xs))
                u(n,m)=rand*(Xs-Xx)+Xx;
            end
        end
    end