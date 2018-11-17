function [A,B] = motionDynamicsLinearize(x,u)
    % A = df/dx(x,u)
    global delx nState delu mControl
    for j = 1:nState
        dx = zeros(nState, 1);
        dx(j) = delx; 
        A(:,j) = (eqnOfMotion2DPointmass(x+dx,u) - eqnOfMotion2DPointmass(x,u))./delx;
    end
    for j = 1:mControl
        du = zeros(mControl, 1);
        du(j) = delu; 
        B(:,j) = (eqnOfMotion2DPointmass(x,u+du) - eqnOfMotion2DPointmass(x,u))./delu;
    end
    
end