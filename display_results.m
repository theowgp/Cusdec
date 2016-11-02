function display_results(solx, solu, N, d, T, x0, v0)


[m, n] = size(solx) ;
mesh = Mesh(T, n-1);
t = mesh.t;

sol = solx';

%% GET ENDTIME VALUES
[xT, vT, zT, uT] = convert(solx(:, end), solu(:, end, 1), N, d);


%% NORM of the SYSTEM VELOCITY at the end-time
normv = norm(solx(N*d+1:2*N*d, end))



%% PLOT THE LYAPUNOV FUNCTION
figure
for k = 1:length(t)
    [x, v, z] = convert_state(sol(k, :), N, d);
    YV(k) =  B(v, N, d);
end
plot(t, YV);
title('B(v, v)');




%% PLOT TRAJECTORIES
figure
for i = 1:N
    plot(sol(1, 2*i-1), sol(1, 2*i), 'o', 'Color', 'k');
    hold all
    plot(sol(:, 2*i-1), sol(:, 2*i), 'Color', 'blue');
    hold all
end
title('evolution');
%% PLOT THE TERMINAL VELOCITY VECTORS
hold all
for i = 1:N
    h = quiver(sol(end, 2*i-1), sol(end, 2*i), 2*vT(i, 1), 2*vT(i, 2),'filled');
    h.MaxHeadSize = 1;
    set(h,'linewidth',1);
    set(h,'color',[1,0,0]);
end
%% PLOT THE INITIAL VELOCITY VECTORS
hold all
for i = 1:N
    h = quiver(x0(i, 1), x0(i, 2), 0.5*v0(i, 1), 0.5*v0(i, 2),'filled');
    h.MaxHeadSize = 0.5;
    set(h,'linewidth',1);
    set(h,'color',[0,0,0]);
end



%% PLOT THE CONTROLS
[mu, nu, su] = size(solu) ;
meshu = Mesh(T, nu-1);
tu = meshu.t;
% d = 1
figure
for i = 1:N
    plot(tu(1:end), solu(2*i-1, :, 1));
    hold all
end
title('controls d = 1');

% d = 2
figure
for i = 1:N
    plot(tu(1:end), solu(2*i, :, 1));
    hold all
end
title('controls d = 2');

% 
% %% PLOT X
% figure
% for k = 1:length(t)
%     x = reshape(sol(k, 1 : N*d), [d, N])';
%     YX(k) =  B(x, x, N);
% end
% plot(t, YX);
% title('X(t)');







BvT = B(vT, N, d)


vT







