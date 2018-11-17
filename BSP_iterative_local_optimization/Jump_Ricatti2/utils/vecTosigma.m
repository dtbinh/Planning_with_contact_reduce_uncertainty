function s = vecTosigma(v, n)
    vtemp = v;
    s1 = zeros(n,n);
    s2 = zeros(n,n);
    for i = 1:n
        s1(i,1:i) = [vtemp(1:i)]';
        s2(i,1:i) = [vtemp(1:i)]';
        s2(i,i) = 0;
        vtemp = vtemp(i+1:end);
    end
    s2 = s2';
    s = s1+s2;
end