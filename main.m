% clear all

%% PARAMETERS:
% number of agents
N = 3;
% dimension
d = 2;

% final time
T = 3;
% number of time windows
ndT = 10;
% number of horizons
nh = 2;

% time window
dT = T/ndT;

% mesh length of a window
n = 10;


% radius of interraction
R = 3;




%% INITIAL CONDITIONS
% initial positions
x0 = initx(N, d, N, R);
% x0 = x00;

% initial velocities
v0 = initv(N, d, 1);
% v0 = v00;


%% CREATE SOLUTION AND CONTROL STACK
solution = zeros(2*N*d+1, (n+1)*ndT);
 
s = 3;% (should be the same as is used in Runge=Kutta scheme)
control = zeros(N*d, n*ndT, s);

% set the initial condition
argx0 = x0;
argv0 = v0;


% set the initial adjacency matrix
Adjc = get_adjacency(x0, N, R, zeros(N));



for k = 1:ndT
    k
    Adjc
    
    [solx, solu, Adjc] = iteration_MPC(N, d, argx0, argv0, dT + nh*dT, n + nh*(n+1), R, Adjc);
    [argx0, argv0, z] = convert_state(solx(:, n+1), N, d);
    
    solution(:, (n+1)*(k-1)+1:(n+1)*k) = solx(:, 1:n+1);
    control(:, n*(k-1)+1:n*k, :) = solu(:, 1:n, :);
    
    plot_trajectories(solx(:, 1:n+1), N);
    hold all
end






[xT, vT] =  display_results(solution, control, N, d, T, x0, v0, Adjc, R);
