function CStep = costFunctionRunning(u,b)
    global Rt Qt nState
    sigma = vecTosigma(b(nState+1:end), nState);
    CStep = u'*Rt*u + trace(sigma*Qt*sigma);
end