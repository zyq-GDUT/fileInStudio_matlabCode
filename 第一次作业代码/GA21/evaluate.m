%��Ӧ������
%fit������ÿ���������Ӧ��
function fit=evaluate(pop,NP)
    fit = zeros(1,NP);
    for i=1:NP    
        fit(i)=evaluationfun(pop(i,:));
    end
end