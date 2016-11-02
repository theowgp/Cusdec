% clear all

%% PARAMETERS:
% number of agents
N = 5;
% dimension
d = 2;

% final time
T = 10;
% number of time windows
ndT = 10;
% time window
dT = T/ndT;

% mesh length of a window
n = 10;


%% SET OBJECTIVE PARAMETERS
alpha1 = 1; % integral of V(t)
alpha2 = 1; % V(T)
alpha3 = 0; % control
alpha5 = 0; % integral of X(t)


%% INITIAL CONDITIONS
% initial positions
x0 = initx(N, d, N);
% x0 = x00;

% initial velocities
v0 = initv(N, d, N);
% v0 = v00;



%% CREATE SOLUTION AND CONTROL STACK
solution = zeros(2*N*d+1, (n+1)*ndT);
 
s = 3;% (should be the same as is used in Runge=Kutta scheme)
control = zeros(N*d, n*ndT, s);

% set the initial condition
argx0 = x0;
argv0 = v0;



for k = 1:ndT
    k
        
    [solx, solu] = Solve(N, d, argx0, argv0, dT, n, alpha1, alpha2, alpha3, alpha5);
    [argx0, argv0, z] = convert_state(solx(:, end), N, d);
    
    solution(:, (n+1)*(k-1)+1:(n+1)*k) = solx;
    control(:, n*(k-1)+1:n*k, :) = solu;
    
    plot_trajectories(solx, N);
    hold all
end








display_results(solution, control, N, d, T, x0, v0);