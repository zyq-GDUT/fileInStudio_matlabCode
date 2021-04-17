function  [] = selection32(NP,u,gen)
%% Ñ¡Ôñ %%
    global x trace Ob fes;
    for m=1:NP
        Ob1(m)=funcl(u(:,m));
    end
    for m=1:NP
        if Ob1(m)<Ob(m)
            x(:,m)=u(:,m);
        end
    end
    for m=1:NP
        Ob(m)=funcl(x(:,m));
        fes=fes+1;
    end
    trace(gen+1)=min(Ob);

