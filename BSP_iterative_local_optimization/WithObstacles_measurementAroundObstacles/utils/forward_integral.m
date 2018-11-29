function [b_new, u_new] = forward_integral(L,I,b_nom,u_nom)
    global nState mControl N x0 sigmanom0 ebs
    
    b_new = zeros(nState+nState*(nState+1)/2, N);
    u_new = zeros(mControl, N);
    
    b_new = zeros(nState+nState*(nState+1)/2, N);
    b_new(:,1) = [x0; sigmanom0];
    for t = 2:N
        u_new(:,t-1) = u_nom(:,t-1) + ebs.*I{t-1,1} + L{t-1,1}*(b_new(:,t-1) - b_nom(:,t-1));
        b_new(:,t) = [eqnOfMotion2DPointmass(b_new(1:nState,t-1), u_new(:,t-1)); motionTerm2(b_new(:,t-1), u_new(:,t-1))];
    end
    
end