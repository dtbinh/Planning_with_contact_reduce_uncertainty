% Cost function for trajectory optimization
function J = costFnTrajOpt(params)
    global N xf nState mControl
    Q = 2*eye(2); R = 1; Qf = Q;
    x = [];
    u = [];
    for i = 1:nState
        x = [x; params(N*(i-1)+1:N*i)];
    end
    for i = 1:mControl
        u = [u; params(nState*N+N*(i-1)+1:nState*N+N*i)];
    end
    Jf = (x(:,end)-xf)'*Qf*(x(:,end)-xf) + u(:,end)'*R*u(:,end);
    Jrun = 0;
    for iter = 1:N-1
        Jrun = Jrun + ((x(:,iter) -xf)'*Q*(x(:,iter) -xf) + u(:,iter)'*R*u(:,iter) );%*dt;  
    end
    J = Jf + Jrun;
end