%% PLOT TRAJECTORIES
for i = 1:N
    plot(sol(:, 2*i-1), sol(:, 2*i));
    hold all
end
title('evolution');

%% PLOT LYAPUNOV FUNCTION
for k = 1:length(t)
    x = reshape(sol(k, 1 : N*d), [d, N])';
    v = reshape(sol(k, N*d+1 : 2*N*d), [d, N])';
    YV(k) =  B(v, v, N);
end

 figure
 plot(t, YV);
 title('V(t)');