function CStep = costFunctionRunning(u,b)
    global Rt Qt nState Ct
    sigma = vecTosigma(b(nState+1:end), nState);
    CStep = u'*Rt*u + trace(sigma*Qt*sigma) + Ct*collisionCost(b(1:nState), sigma);
end