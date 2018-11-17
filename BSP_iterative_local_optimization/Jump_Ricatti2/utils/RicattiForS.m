function [S1, s2, s3, L, I] = RicattiForS()
    global F G Fi Gi ei Qt_full P R q r p nState mControl N S1T s2T s3T
    
%     S1 = zeros(N,nState+nState*(1+nState)/2, nState+nState*(1+nState)/2);
%     s2 = zeros(N,nState+nState*(1+nState)/2, 1);
%     s3 = zeros(N,1);
%     L = zeros(N, mControl, nState+nState*(1+nState)/2);
%     I = zeros(N,mControl,1);
    
    S1 = cell(N,1);
    s2 = cell(N,1);
    s3 = cell(N,1);
    L = cell(N,1);
    I = cell(N,1);
    
    S1{N,1} = S1T;
    s2{N,1} = s2T;
    s3{N,1} = s3T;

    for t = N-1:-1:1
        Cterm = zeros(nState+nState*(1+nState)/2,nState+nState*(1+nState)/2);
        cterm = zeros(nState+nState*(1+nState)/2,1);
        Dterm = zeros(mControl,mControl);
        dterm = zeros(mControl,1);
        Eterm = zeros(mControl, nState+nState*(1+nState)/2);
        eterm = 0;
        for i = 1:nState
%             Fitrans = permute(Fi(t,:,:,i), [1 3 2 4]);
%             Fitrans = squeeze(Fitrans);
%             Gitrans = permute(Gi(t,:,:,i), [1 3 2 4]);
%             Gitrans = squeeze(Gitrans);
            Cterm = Cterm + Fi{t,i}'*S1{t+1,1}*Fi{t,i};
            cterm = cterm + Fi{t,i}'*S1{t+1,1}*ei{t,i};
            Dterm = Dterm + Gi{t,i}'*S1{t+1,1}*Gi{t,i};
            dterm = dterm + Gi{t,i}'*S1{t+1,1}*ei{t,i};
            Eterm = Eterm + Gi{t,i}'*S1{t+1,1}*Fi{t,i};
            eterm = eterm + 0.5*ei{t,i}'*S1{t+1,1}*ei{t,i};
        end
        C = Qt_full{t,1} + F{t,1}'*S1{t+1,1}*F{t,1} + Cterm;
        c = q{t,1} + F{t,1}'*s2{t+1,1} + cterm;
        D = R{t,1} + G{t,1}'*S1{t+1,1}*G{t,1} + Dterm;
        d = r{t,1} + G{t,1}'*s2{t+1,1} + dterm;
        E = P{t,1} +  G{t,1}'*S1{t+1,1}*F{t,1} + Eterm;
        e = p{t,1} + s3{t+1,1} + eterm;
        
        S1{t,1} = C - E'*inv(D)*E;
        s2{t,1} = c - E'*inv(D)*d;
        s3{t,1} = e - 0.5*d'*inv(D)*d;
        L{t,1} = D\E.*(-1);
        I{t,1} = D\d.*(-1);
    end
end