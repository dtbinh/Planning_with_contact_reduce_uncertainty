% EKF
clc;
clear all;
close all;

global x0 dt k m w x_traj t_traj nState mControl pMeasure R Q C D

x0 = 0.1;
T = 10; %sec
N = 1000;
dt = T/N;
nState = 2;
mControl = 1;
pMeasure = 2;

% Constants
k = 100;
m = 10;
w = sqrt(k/m);

% Motion
R = x0*ones(nState);

% Measurement
C = eye(pMeasure);
D = zeros(pMeasure, mControl);
Q = x0*eye(nState);

% Generating actual trajectory in advance for 10s using ODE45
xinit = [x0;0];
t0 = linspace(0,T,N);
[t_traj,x_traj] = ode45(@DynamicsODE45, t0, xinit);
x_traj = x_traj';

%% Initial estimate
mu(:,1) = [0;0];
sigma(:,:,1) = x0.*zeros(nState, nState);
count = 1;
for t = 0:dt:T-dt
    
    mu_bar = eulerIntegrationMotion(mu(:,count), 0);
    [G, H] = linearize();
    sigma_bar = G*sigma(:,:,count)*G' + R;
    K = sigma_bar*H'*inv(H*sigma_bar*H' + Q);
    z = dynamics_spring_mass(t+dt);
    mu(:, count+1) = mu_bar + K*(z - measurementDynamics(mu_bar,0));
    sigma(:,:,count +1) = (eye(nState) - K*H)*sigma_bar;
    count = count + 1;
end

figure;

plot(t0, x_traj(1,1:1000));
hold on;
plot(t0, mu(1,1:1000));

% Linearize dynamics and measurement for a spring mass system
function [G, H] = linearize()
    global nState k m mControl dt C D pMeasure
    A = [0 1;-k/m 0];
    B = [0;0];
    M = [A B; zeros(1,nState+mControl)].*dt;
    MM = expm(M);
    G = MM(1:nState,1:nState);
    Mz = [C D; zeros(1,nState+mControl)].*dt;
    MMz = expm(Mz);
    H = MMz(1:pMeasure,1:nState);
end

% Euler integral for one step motion
function mut = eulerIntegrationMotion(x, u)
    global k m dt nState
    mut = zeros(nState, 1);
    mut(1) = x(1) + x(2)*dt;
    mut(2) = x(2) +(-k/m*x(1)+u)*dt;
end

% Mapping from current belief to the measurement
function z = measurementDynamics(mu_bar,u)
    global C D pMeasure
    z = zeros(pMeasure, 1);
    z = C*mu_bar + D*u;
end

% Deterministic measurement
function z = dynamics_spring_mass(t)
    global x0 w pMeasure
    z = zeros(pMeasure, 1);
    % Equation of motion
    z(1) = x0*cos(w*t);
    z(2) = -x0*w*sin(w*t);
    
end

function p_zx = observation_Likelihood(measurement, X)
    mux = measurement(1);
    sigmax = 0.005;
    muy = measurement(2);
    sigmay = 0.1;
%     p_zx = 1/(sigmax*sqrt(2*pi))*exp( -(((X(1)-mux)^2)/(2*sigmax^2)) );
    p_zx = 1/(sigmax*sigmay*2*pi)*exp( -(((X(1)-mux)^2)/(2*sigmax^2) + ((X(2)-muy)^2)/(2*sigmay^2)) );
end

function dxdt = DynamicsODE45(t,x)
    global k m
    dxdt(1) = x(2);
    dxdt(2) = -k/m*x(1);
    dxdt = dxdt';
end