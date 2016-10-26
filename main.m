
%% PARAMETERS:
% number of agents
N = 3;
% dimension
d = 2;
% final time
T = 100;
% mesh length
n = 800;
% create mesch
mesh = Mesh(T, n);



%% INITIAL CONDITIONS
% initial positions
x0 = initx(N, d, N);

% initial velocities
v0 = initv(N, d, 2);





%% CREATE THE DYNAMICS
delta = 1;
dynamics = Dynamics(N, d, delta);



%% CREATE THE OBJECTIVE
objective = Objective(N, d);



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


%% MINIMIZATION
eps = 1;% not used 
sigma = 0.001;
limitLS = 5;
limitA = 10;
[solx, solu] = NCG(rk, objective, mesh, solu0, eps, sigma, limitLS, limitA);

sol = solx';
t = mesh.t;





 
 
OUTPUT % script