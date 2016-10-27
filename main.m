% clear all

%% PARAMETERS:
% number of agents
N = 3;
% dimension
d = 2;

% final time
T = 20;
% number of time windows
ndT = 20;
% time window
dT = T/ndT;

% mesh length of a window
n = 40;


% radius of interraction
R = 3;



%% INITIAL CONDITIONS
% initial positions
initial_x0 = initx(N, d, N);
% x0 = x00;

% initial velocities
initial_v0 = initv(N, d, N);
% v0 = v00;


%% CREATE SOLUTION AND CONTROL STACK
solution = zeros(2*N*d+1, (n+1)*ndT);
 
s = 3;% (should be the same as is used in Runge=Kutta scheme)
control = zeros(N*d, n*ndT, s);

% set the initial condition
argx0 = initial_x0;
argv0 = initial_v0;

for k = 1:ndT
    [solx, solu] = iteration_MPC(N, d, argx0, argv0, dT, n, R);
    [argx0, argv0, z] = convert_state(solx(:, end), N, d);
    
    solution(:, (n+1)*(k-1)+1:(n+1)*k) = solx;
    control(:, n*(k-1)+1:n*k, :) = solu;
    
    plot_trajectories(solx, N);
    hold all
end

% sol = solution';
% t = mesh.t;




% OUTPUT; % script