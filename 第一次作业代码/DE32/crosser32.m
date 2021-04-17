function u = crosser32(D,CR,v,x,NP,Xx,Xs)
%% ½»²æ %%
    r=randi([1,D],1,1);
    for n=1:D
        cr=rand(1);
        if(cr<=CR)||(n==r)
            u(n,:)=v(n,:);
        else
            u(n,:)=x(n,:);
        end
    end
    %% ±ß½çÎüÊÕ %%
    for n=1:D
        for m=1:NP
            if u(n,m)<Xx
                u(n,m)=Xx;
            end
            if u(n,m)>Xs
                u(n,m)=Xs;
            end
        end
    end

