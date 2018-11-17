function CostT = costTerminal(b)
    global Rt Qf nState xf
    x = b(1:nState);
    
    sigma = vecTosigma(b(nState+1:end), nState);
    CostT = (x-xf)'*Qf*(x-xf)+ trace(sigma*Qf*sigma);
end