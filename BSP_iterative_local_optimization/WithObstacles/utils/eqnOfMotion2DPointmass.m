%Eqn of motion for 2D point mass system x(t+1) = f(xt, ut)
function f = eqnOfMotion2DPointmass(x, u)
global dt
    f = x + u*dt; % u is the velocity
end