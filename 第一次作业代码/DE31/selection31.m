function  [] = selection31(NP,u,gen)
%% ѡ?????? %%
    global x Ob trace31 fes;
    Ob1 = zeros(NP);
    for m=1:NP
        Ob1(m)=funcl(u(:,m));
    end
    for m=1:NP
        if Ob1(m)<Ob(m)
            x(:,m)=u(:,m);
        end
    end
    for m =1:NP
        Ob(m)=funcl(x(:,m));
        fes=fes+1;
    end
    trace31(gen+1)=min(Ob(1));
    

