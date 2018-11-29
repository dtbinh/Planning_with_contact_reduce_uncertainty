function s30 = candidateCost(L, b_new, u_new)
    global nState mControl N
    
    S1 = cell(N,1);
    s2 = cell(N,1);
    s3 = cell(N,1);
    
    %%Finding terminal S (constant)
    [S1T, s2T, s3T] = terminalS(b_new, u_new);
    
    S1{N,1} = S1T;
    s2{N,1} = s2T;
    s3{N,1} = s3T;
    
    % Linearize motion dynamics
    [F, G] = stochasticMotionDynamicsLinearize(b_new, u_new);
    % Linearize observation dynamics
    [Fi, Gi, ei] = stochaticVarianceLinearize(b_new, u_new);
    % Quadratize cost function
    [Qt_full, P, R, q, r, p] = quadratizeCost(b_new, u_new);
    
    for t = N-1:-1:1
        Sterm = zeros(nState+nState*(1+nState)/2,nState+nState*(1+nState)/2);
        sterm = 0;
        for i = 1:nState
            Sterm = Sterm + (Fi{t,i} + Gi{t,i}*L{t,1})'*S1{t+1,1}*(Fi{t,i} + Gi{t,i}*L{t,1});
            sterm = sterm + 0.5*ei{t,i}'*S1{t+1,1}*ei{t,i};
        end
        S1{t,1} = Qt_full{t,1} + L{t,1}'*R{t,1}*L{t,1} + L{t,1}'*P{t,1} + P{t,1}'*L{t,1} + (F{t,1} + G{t,1}*L{t,1})'*S1{t+1,1}*(F{t,1} + G{t,1}*L{t,1}) + Sterm;
        s3{t,1} = p{t,1} + s3{t+1,1} + sterm;
    end
    
    s30 = s3{1,1};
end