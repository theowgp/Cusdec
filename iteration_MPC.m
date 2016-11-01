function [solx, solu, Adjc] = iteration_MPC(N, d, argx0, argv0, T, n, R, Adjc)

    Adjc =  get_adjacency(argx0, N, R, Adjc);

                     
    [solx, solu] = Solve(N, d, argx0, argv0, T, n, R, Adjc);

end






