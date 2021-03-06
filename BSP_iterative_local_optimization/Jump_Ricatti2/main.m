% 2D point robot moving in an obstacle free space (velocity is the control input, no inertia system)
% assuming motion noise (Rt) to be independent of x and u (unit variance)
% measurement noise (Qt) is dependent on x as given by the contour plot
% variance is stored rowwise and its square-root is taken elementwise
% obstacles are square and 2D, stored in anticlockwise order starting from bottom-left

clc;
clear all;
close all;
addpath('/media/saumya/Data/Study/CMU/Robotics/PlanningWithContact/Codes/Planning_with_contact_to_reduce_uncertainty/BSP_iterative_local_optimization/Jump_Ricatti2/utils');
global nState mControl pMeasure xf x0 N dt umax T Amp sig mu F G Fi Gi ei Qt_full P R q r p mapxmax mapymax delx delu delsig Mt Rt Qf Qt S1T s2T s3T x_nom0 sigmanom0 u_nom0 obstacles ebs

%% Declaring parameters
% Q = 2*eye(2); R = 1; Qf = Q;
nState = 2; % size of state space
mControl = 2; % size of control space
pMeasure = 2;

N = 100 ;%number of grid points (time steps)
dt = 0.1 ;%time step for Euler Integration
T = N*dt ;%total time
t0 = [linspace(0,T,N)]';

x0 = [10;20];
xf = [20;70];
umax = inf ;%upper control bound (optimization variable)

% Parameters for the onbservation map
Amp = 0.002;
sig = 10000;
mu = 60;

% Dimensions of the map
mapxmax = 100;
mapymax = 100;

% delta steps for finite difference operations
delx = 0.01;
delu = 0.01;
delsig = 0.01;

% Noise for the motion model
Mt = 1*eye(nState);

% Cost funcition penalty terms
Rt = 1*eye(mControl); %u'*Rt*u
Qt = 400*eye(nState);  % sigma*Qt*sigma
Qf = 100*N*eye(nState); %x'*Qf*x + sigma*Qf*sigma

% Obstacles are 
obstacles{1} = [30 20; 60 20; 60 80; 30 80]';

% Increment factor for the control
ebs = 1;

%% Loading nominal trajectory
load nominalTraj_X_origin15_15.mat
load nominalTraj_U_origin15_15.mat
load nominalTraj_t0_origin15_15.mat

%% Setting up fmincon to find nominal trajectory
% A = [] ;%empty because no linear equations
% b = [] ;%empty because no linear equations
% Aeq = [] ;%empty because no linear equations
% beq = [] ;%empty because no linear equations
% options =optimoptions(@fmincon,'TolFun',0.0001,'MaxIter',10000,...
%     'MaxFunEvals',100000,'Display','iter','DiffMinChange',0.001,'Algorithm','sqp');
% % State and control bounds
% UB = [inf*ones(1,nState*N),umax*ones(1,mControl*N)];
% LB = -UB;
% params0 = [zeros(1,nState*N), zeros(1,mControl*N)];
% % Finding the nominal trajectory
% COSTFUN = @(params) costFnTrajOpt(params);
% CONSTRAINTFUN = @(params) nonlinconstraintsDirCol(params,@Dynamics_2DpointBot);
% params = fmincon(COSTFUN,params0,A,b,Aeq,beq,LB,UB,CONSTRAINTFUN,options);
% x = [];
% u = [];
% for i = 1:nState
%     x = [x; params(N*(i-1)+1:N*i)];
% end
% for i = 1:mControl
%     u = [u; params(nState*N+N*(i-1)+1:nState*N+N*i)];
% end
% save('nominalTraj_X_origin15_15.mat', 'x');
% save('nominalTraj_U_origin15_15.mat', 'u');
% save('nominalTraj_t0_origin15_15.mat', 't0');
x_nom0 = x;
u_nom0 = u;

%% Plotting initial trajectory and environment
figure(1);
set(gcf,'units','points','position',[10,10,700,700])
% Observation environment variance
subplot(4,2,[5,6,7,8])
mapx = 0:mapxmax;
mapy = 0:mapymax;
[mapX,mapY] = meshgrid(mapx,mapy);

% contZvar = Amp*exp(-(mapX - mu).^2/sig);
contZvar = Amp*(mu - mapX).^2;
conto = contourf(mapX,mapY,contZvar,50, 'edgecolor','none');
colormap(flipud(gray));
hold on;
% Plotting nominal trajectory
% plotNominalTraj(t0, x_nom0, u_nom0);
hold on;

%% Initializing Nominal belief trajectory
sigmanom0 = sigmaToVec(50*eye(nState), nState);
b_nom = zeros(nState+nState*(nState+1)/2, size(x,2)); % zero variance
b_nom(1:nState, :) = x_nom0;
b_nom(nState+1:end, 1) = sigmanom0;
u_nom = u_nom0;


outer_iter = 1;
max_iterations_outer = 20;

% Create a movie structure to add frames to
mov(1:max_iterations_outer) = struct('cdata', [], 'colormap', []);
% Create a video writer to use the writeVideo function
v = VideoWriter('BSP_iterative_local_optimization_MeanJumpRicatti.avi');
% Make the video writer available for writing
open(v);
TotCostTraj = 100000000;
%% iLQG loop
while(outer_iter <= max_iterations_outer) % or l(t) <lt
    if ebs == 1 % Plot only if new trajectory is accepted
        plotTrajectoryWithVariance(t0, b_nom, u_nom);
    end
    %%Finding terminal S (constant)
    [S1T, s2T, s3T] = terminalS(b_nom, u_nom);
    
    % Linearize motion dynamics
    [F, G] = stochasticMotionDynamicsLinearize(b_nom, u_nom);
    % Linearize observation dynamics
    [Fi, Gi, ei] = stochaticVarianceLinearize(b_nom, u_nom);
    % Quadratize cost function
    [Qt_full, P, R, q, r, p] = quadratizeCost(b_nom, u_nom);
    % Ricatti for S
    [S1, s2, s3, L, I] = RicattiForS();
    
    % ODE45 for forward propagation
    [b_new, u_new] = forward_integral(L,I,b_nom,u_nom);
    
%     %% WRONG
%     TotCostTraj_new = 0;
%     for i = 1:N-1
%         TotCostTraj_new = TotCostTraj_new + costFunctionRunning(u_new(:,i),b_new(:,i));
%     end
%     TotCostTraj_new = TotCostTraj_new + costTerminal(b_new(:,N));
%     %%
%     if TotCostTraj_new < TotCostTraj
%         b_nom = b_new;
%         u_nom = u_new;
%         TotCostTraj = TotCostTraj_new;
%         ebs = 1;
%     else
%         ebs = ebs/2;
%     end

     b_nom = b_new;
     u_nom = u_new;
    
    if ebs ==1 % Store frame only if new trajectory is accepted
        % Video writing
        ax = gcf();
        mov(outer_iter) = getframe(ax);
        % Write frame to the video writer "v"
        writeVideo(v,mov(outer_iter));

        outer_iter = outer_iter +1;
    end
    pause(0.01);

end
close(v);