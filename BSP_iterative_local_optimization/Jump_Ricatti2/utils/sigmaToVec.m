function vec = sigmaToVec(sig, n)
%     vec = zeros(n*(n+1)/2,1);
    vec = [];
    for i = 1:n
        vec = [vec; sig(i, 1:i)'];
    end
end