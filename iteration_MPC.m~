function [solx, solu, Adjc] = iteration_MPC(N, d, argx0, argv0, T, n, R, Adjc)


    %% SET ADJACENCY MARIX
    Adjc =  get_adjacency(argx0, N, R, Adjc);
%     Adjc = ones(N) - eye(N);
    

    %% SET ARRAY OF NUMBERS OF AGENTS FOR EACH OPC
    Ns = zeros(1, N);

    %% SET SOLUTION VECTOR
    solx = zeros(2*N*d+1, n+1);

    %% SET THE CONTROL VECTOR 
    s = 3; %(should be the same as is used in Runge=Kutta scheme)
    solu = zeros(N*d, n, s);
% 
for p = 1:N
%     p
       
    x0 = set_w0(argx0, Adjc, N, p);
    v0 = set_w0(argv0, Adjc, N, p);
    Ns(p) = numbers_N(Adjc, N, p); 
    
       
                     
    [solxp, solup] = Solve(Ns(p), d, x0, v0, T, n, R);

%     stack on the positions
    solx((p-1)*d+1:p*d, :) = solxp(1:d, :);

%     stack on the velocities
    solx(N*d+(p-1)*d+1:N*d+p*d, :) = solxp(Ns(p)*d+1:Ns(p)*d+d, :);

%     stack on the controls
    solu((p-1)*d+1:p*d, :, :) = solup;

end





end

