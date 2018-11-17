function [F, G] = stochasticMotionDynamicsLinearize(b_nom, u_nom)
    global nState mControl delx delu delsig N
    % F = [F1 F2; F3 F4]
    Ftemp = zeros(nState+nState*(1+nState)/2,nState+nState*(1+nState)/2);
    F  = cell(N,1);
    F1 = zeros(nState, nState);
    F2 = zeros(nState, nState*(1+nState)/2);
    F3 = zeros(nState*(1+nState)/2, nState);
    F4 = zeros(nState*(1+nState)/2, nState*(1+nState)/2);
    
    Gtemp = zeros(nState+nState*(1+nState)/2,mControl);
    G = cell(N,1);
    
    for i = 1:N
        x = b_nom(1:nState, i);
        sigma = b_nom(nState+1:end, i);
        u = u_nom(:,i);
        %Finding F
        F1 = motionDynamicsLinearize(x,u);
        for j = 1:nState
            dx = zeros(nState, 1);
            dx(j) = delx;
            F3(:,j) = (motionTerm2([x+dx;sigma],u) - motionTerm2([x;sigma],u))./delx;
        end
        
        for j = 1:nState*(nState+1)/2
            dsig = zeros(nState*(nState+1)/2, 1);
            dsig(j) = delsig;
            F4(:,j) = (motionTerm2([x; sigma+dsig],u) - motionTerm2([x;sigma],u))./delsig;
        end
        
        Ftemp(1:nState, 1:nState) = F1;
        Ftemp(1:nState, nState+1:end) = F2;
        Ftemp(nState+1:end, 1:nState) = F3;
        Ftemp(nState+1:end, nState+1:end) = F4;
        F{i} = Ftemp;
        
        % Finding G
        [~,B] = motionDynamicsLinearize(x,u);
        Gtemp(1:nState,1:mControl) = B;
        for j = 1:mControl
            du = zeros(mControl, 1);
            du(j) = delu; 
            G_(:,j) = (motionTerm2([x;sigma],u+du) - motionTerm2([x;sigma],u))./delu;
        end
        Gtemp(nState+1:end,:) = G_;
        G{i} = Gtemp;
    end
    
end