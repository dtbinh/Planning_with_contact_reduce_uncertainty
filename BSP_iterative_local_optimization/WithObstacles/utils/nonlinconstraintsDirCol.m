% Non-linear constraints for direct collocation
function [c,ceq] = nonlinconstraintsDirCol(params,Dynamics) 
    global N x0 xf dt nState mControl
    x = [];
    u = [];
    for i = 1:nState
        x = [x; params(N*(i-1)+1:N*i)];
    end

    for i = 1:mControl
        u = [u; params(nState*N+N*(i-1)+1:nState*N+N*i)];
    end
    
    ceq_0 = x(:,1)-x0;
    ceq_f = x(:,end)-xf;
    ceq = ceq_0;
    x_dot_k = Dynamics(x(:,1),u(:,1));
    
    for k = 1:N-1        
        x_k = (x(:,k));
        x_k_p1 = x(:,k+1);
        x_dot_k_p1 = Dynamics(x(:,k+1),u(:,k+1));

        h = dt;
        x_ck = 1/2*(x_k + x_k_p1) + h/8*(x_dot_k - x_dot_k_p1);
        u_ck = (u(:,k)+u(:,k+1))/2;
        x_dot_ck = Dynamics(x_ck,u_ck);
      
        defects = (x_k - x_k_p1) + dt/6* (x_dot_k + 4*x_dot_ck + x_dot_k_p1);
        ceq = [ceq;defects];
        x_dot_k = x_dot_k_p1;
    end
    c = [] ;%inequality <= constraint
    ceq = [ceq;ceq_f] ;%equality = constraint
end