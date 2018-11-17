function W = stochasticVarianceTermW(b,u)
    global nState Mt
    [A,~] = motionDynamicsLinearize(b(1:nState),u);
    H = findingHLinearizeMeasurement(b(1:nState),u);
    N = findingN_MeasurementVariance(b(1:nState),u);
    
    sigma = vecTosigma(b(nState+1:end), nState);
    tau  = A*sigma*(A*sigma)' + Mt*Mt';
    K = tau*H'*inv(H*tau*H' + N*N');
    
    W = (K*H*tau).^0.5;
    W = [W; zeros(nState*(nState+1)/2, nState)];
end