function result=func1(x) %result需要接收返回值
    summ=sum(x.^2);%默认是按列每个元素的平方求和，sum(x.^2,2)则是按行累加
    result=summ;
end