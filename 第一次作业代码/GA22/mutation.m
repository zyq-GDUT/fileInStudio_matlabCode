function mpop = mutation(cpop, Pm)
    global NP D Xx Xs
    for m=1:NP
        for n=1:D
            r=rand(1,1);
            if r < Pm
                cpop(n,m) = rand(1,1)*(Xs-Xx)+Xx;
            end
        end
    end
    mpop = cpop;
end

