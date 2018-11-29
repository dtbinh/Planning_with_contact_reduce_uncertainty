function mes = measurementVarianceMap(x,y)
    global obstacles
    slope =0.3;
    mes = 4 + (1/(1+exp(slope*(x-obstacles{1}(1,1)))) + 1/(1+exp(-slope*(x-obstacles{1}(1,2))))-1)*(1/(1+exp(-slope*(y-obstacles{1}(2,1)))) + 1/(1+exp(slope*(y-obstacles{1}(2,3))))-1);
end