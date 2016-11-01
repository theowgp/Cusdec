function [solxp, solup] = Solve(N, d, x0, v0, T, n, R, Adjc)




%% SET OBJECTIVE PARAMETERS
alpha1 = 0; % integral of V(t)
alpha2 = 1; % V(T)
alpha3 = 0; % control
alpha5 = 0; % integral of X(t)
alpha7 = 1; % the potential 

gamma = 0;
M = 1;



%% CREATE THE DYNAMICS
delta = 1;
dynamics = Dynamics(N, d, gamma, delta, alpha1, alpha3, alpha5, alpha7, M, R, Adjc);

%% CREATE THE OBJECTIVE
objective = Objective(dynamics, N, d,  alpha2);



%% CREATE RUNGE KUTTA SOLVER
A = [0 0 0; 0.5 0 0; -1 2 0];
b = [1.0/6.0    2.0/3.0    1.0/6.0];
% c = [0  0.5  1];
s = 3;

rkNu = N*d;

arg0 = [reshape(x0', [N*d, 1]); reshape(v0', [N*d, 1]); 0];

rkN = 2*N*d+1;

rk = RungeKutta(A, b, s, dynamics, objective, arg0, rkN, rkNu, T, n);


%% INITIAL CONTROL GUESS
solu0 = zeros(rkNu, n,  s);


%% NCG MINIMIZATION
eps = 1;% not used 
sigma = 0.001;
limitLS = 100;
limitA = 25;

mesh = Mesh(T, n);

[solxp, solup] = NCG(rk, objective, mesh, solu0, eps, sigma, limitLS, limitA);


end

