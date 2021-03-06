function [solx, solu, Adjc] = iteration_MPC(N, d, argx0, argv0, T, n, R, Rh, Adjc, decentralized)


    %% SET ADJACENCY MARIX
    Adjc =  get_adjacency(argx0, N, R, Rh, Adjc);
%     Adjc = ones(N) - eye(N);
    

    %% SET ARRAY OF NUMBERS OF AGENTS FOR EACH OPC
    Ns = zeros(1, N);

    %% SET SOLUTION VECTOR
    solx = zeros(2*N*d+1, n+1);

    %% SET THE CONTROL VECTOR (should be the same as is used in Runge=Kutta scheme)
    s = 3;
    solu = zeros(N*d, n, s);

if decentralized
    for i = 1:N
%         i
       
        x0 = set_w0(argx0, Adjc, N, i);
        v0 = set_w0(argv0, Adjc, N, i);
        Ns(i) = numbers_N(Adjc, N, i); 
        
        Adjci = zeros(Ns(i));
        Adjci(1, 2:Ns(i)) = ones(1, Ns(i)- 1);
        Adjci(2:Ns(i), 1) = ones(Ns(i)- 1, 1);
        
        
        [solxi, solui] = Solve(Ns(i), d, x0, v0, T, n, R, Adjci);

    %       for the first dimension
        solx((i-1)*d+1:i*d, :) = solxi(1:d, :);
    %       for the second dimension
        solx(N*d+(i-1)*d+1:N*d+i*d, :) = solxi(Ns(i)*d+1:Ns(i)*d+d, :);
    %       for the first dimension
        solu((i-1)*d+1:i*d, :, :) = solui(1:d, :, :);

    end
else
    [solx, solu] = Solve(N, d, argx0, argv0, T, n, R, Adjc);
end





end

