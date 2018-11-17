function H = findingHLinearizeMeasurement(x,u)
    % finite differential of the measurement function
    global nState pMeasure delx
    H = zeros(pMeasure, nState);
    for i = 1:nState
        dx = zeros(nState, 1);
        dx(i) = delx; 
        H(:,i) = (measurement(eqnOfMotion2DPointmass(x+dx,u)) - measurement(eqnOfMotion2DPointmass(x,u)))./delx;
    end
end