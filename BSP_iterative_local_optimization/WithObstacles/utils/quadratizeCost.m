function [Qt_full, P, R, q, r, p] = quadratizeCost(b_nom, u_nom)
    global N nState mControl Qt Rt delx delu delsig
%     q = zeros(N, nState+nState*(1+nState)/2,1);
%     r = zeros(N, mControl,1);
%     p = zeros(N,1);
%     Qt_full = zeros(N, nState+nState*(1+nState)/2, nState+nState*(1+nState)/2);
%     R = zeros(N, mControl,mControl);
%     P = zeros(N, mControl, nState+nState*(1+nState)/2);
    
    q = cell(N,1);
    r = cell(N,1);
    Qt_full = cell(N,1);
    R = cell(N,1);
    P = cell(N,1);
    p = cell(N,1);
    
    Q_ = zeros(nState*(1+nState)/2, nState*(1+nState)/2);

    for t = 1:N
        b = b_nom(:, t);
        u = u_nom(:,t);
        Qt_full{t,1} = zeros(nState+nState*(1+nState)/2, nState+nState*(1+nState)/2);
        for i = 1:nState*(1+nState)/2
            for j = 1:nState*(1+nState)/2
                dbi = zeros(nState+nState*(1+nState)/2,1);
                dbj = zeros(nState+nState*(1+nState)/2,1);
                dbi(nState + i) = delsig;
                dbj(nState + j) = delsig;
                Q_(i,j) = (costFunctionRunning(u,b+dbi+dbj) - costFunctionRunning(u,b-dbi+dbj) - costFunctionRunning(u,b+dbi-dbj) + costFunctionRunning(u,b-dbi-dbj))/(4*delsig*delsig);
            end
        end
        
        Qt_full{t}(nState+1:end, nState+1:end) = Q_;
        
        R{t} = Rt + Rt';
        for i = 1:mControl
            du = zeros(mControl, 1);
            du(i) = delu;
            r{t,1}(i,1) = (costFunctionRunning(u+du,b) - costFunctionRunning(u,b))/delu;
            for j = 1:nState+nState*(1+nState)/2
                db = zeros(nState+nState*(1+nState)/2, 1);
                db(j) = delx;
                % Central second partial derivative
                P{t,1}(i,j) = (costFunctionRunning(u+du,b+db) - costFunctionRunning(u+du,b-db) - costFunctionRunning(u-du,b+db) + costFunctionRunning(u-du,b-db))/(4*delu*delx);
            end
        end
        
        for j = 1:nState+nState*(1+nState)/2
            db = zeros(nState+nState*(1+nState)/2, 1);
            db(j) = delx;
            % Central second partial derivative
            q{t,1}(j,1) = (costFunctionRunning(u,b+db) - costFunctionRunning(u,b))/delx;
        end
        p{t,1} = costFunctionRunning(u,b);
    end
    

end