function N = findingN_MeasurementVariance(x,u)
    global Amp mu sig
    xtp1 = eqnOfMotion2DPointmass(x, u);
%     contZvar = Amp*exp(-(xtp1(1) - mu).^2/sig);
    contZvar = Amp*(mu - xtp1(1)).^2;
    N = [contZvar 0; 0 contZvar];
end