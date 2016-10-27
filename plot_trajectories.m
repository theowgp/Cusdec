function  plot_trajectories(solx, N)
sol = solx';

%% PLOT TRAJECTORIES
for i = 1:N
    plot(sol(:, 2*i-1), sol(:, 2*i));
    hold all
end



end

