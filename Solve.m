function [solx, solu] = Solve(N, d, x0, v0, T, n, alpha1, alpha2, alpha3, alpha5)





%% CREATE THE DYNAMICS
delta = 1;
dynamics = Dynamics(N, d, delta, alpha1, alpha3, alpha5);


%% CREATE THE OBJECTIVE
objective = Objective(dynamics, N, d, alpha2);



%% CREATE RUNGE KUTTA SOLVER
A = [0 0 0; 0.5 0 0; -1 2 0];
b = [1.0/6.0    2.0/3.0    1.0/6.0];
% c = [0  0.5  1];
s = 3;
Nu = N*d;

arg0 = [reshape(x0', [N*d, 1]); reshape(v0', [N*d, 1]); 0];

rk = RungeKutta(A, b, s, dynamics, objective, arg0, 2*N*d+1, Nu, T, n);


%% INITIAL CONTROL GUESS
solu0 = zeros(N*d, n,  s);


%% NCG MINIMIZATION
eps = 1;% not used 
sigma = 0.001;
limitLS = 15;
limitA = 25;
mesh = Mesh(T, n);
[solx, solu] = NCG(rk, objective, mesh, solu0, eps, sigma, limitLS, limitA);



end

