function [S1T, s2T, s3T] = terminalS(b_nom, u_nom)
    global nState Qf delx delsig
    
    b = b_nom(:,end);
    u = u_nom(:,end);
    S1T = zeros(nState+nState*(1+nState)/2, nState+nState*(1+nState)/2);
    s2T = zeros(nState+nState*(1+nState)/2,1);
    s3T = 0;
    
    S11 = (Qf+Qf');
    S21 = zeros(nState*(1+nState)/2, nState);
    S12 = S21';
    S22 = zeros(nState*(1+nState)/2,nState*(1+nState)/2);
    
    for i = 1:nState*(1+nState)/2
        for j = 1:nState*(1+nState)/2
            dbi = zeros(nState+nState*(1+nState)/2,1);
            dbj = zeros(nState+nState*(1+nState)/2,1);
            dbi(nState + i) = delsig;
            dbj(nState + j) = delsig;
            S22(i,j) = (costTerminal(b+dbi+dbj) - costTerminal(b-dbi+dbj) - costTerminal(b+dbi-dbj) +costTerminal(b-dbi-dbj))/(4*delsig*delsig);
%             S22(i,j) = (costTerminal(b+dbi+dbj) - costTerminal(b+dbi) - costTerminal(b+dbj) +costTerminal(b))/(delsig*delsig);

        end
    end
    
    S1T = [S11 S12; S21 S22];
    

    for j = 1:nState+nState*(1+nState)/2
        db = zeros(nState+nState*(1+nState)/2, 1);
        db(j) = delx;
        % Central second partial derivative
        s2T(j,1) = (costTerminal(b+db) - costTerminal(b))/delx;
    end
    s3T = costTerminal(b);
end