function v = mutation31(gen,G,F0,x,NP)
%% £®◊‘  ”¶£©±‰“Ï %%
    lamda = exp(1-G/(G+1-gen));
    F= F0*2^(lamda);
    for m=1:NP
        r1=randi([1,NP],1,1);
        while (r1==m) 
            r1=randi([1,NP],1,1);
        end
        r2=randi([1,NP],1,1);
        while(r2==m)||(r2==r1)
            r2=randi([1,NP],1,1);
        end
        r3=randi([1,NP],1,1);
        while(r3==m)||(r3==r1)||(r3==r2)
            r3=randi([1,NP],1,1);
        end
        v(:,m)=x(:,r1)+F*(x(:,r2)-x(:,r3));
    end

