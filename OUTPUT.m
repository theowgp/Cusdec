mesh = Mesh(dT, n);

sol = solution';
t = zeros(1, ndT*(n+1));
for k = 1:ndT 
    t((n+1)*(k-1)+1:(n+1)*k) = mesh.t + (k-1)*dT;
end


%% GET ENDTIME VALUES
[xT, vT, zT, uT] = convert(solx(:, end), solu(:, end, 1), N, d);


%% NORM of the SYSTEM VELOCITY at the end-time
normv = norm(solx(N*d+1:2*N*d, end))



%% PLOT THE LYAPUNOV FUNCTION
figure
for k = 1:length(t)
%     x = reshape(sol(k, 1 : N*d), [d, N])';
    v = reshape(sol(k, N*d+1 : 2*N*d), [d, N])';
    YV(k) =  B(v, v, N);
end
plot(t, YV);
title('V(t) = 1/2N^2  sumij||vi -vj ||^2');




%% PLOT TRAJECTORIES
figure
for i = 1:N
    plot(sol(1, 2*i-1), sol(1, 2*i), 'o', 'Color', 'k');
    hold all
    plot(sol(:, 2*i-1), sol(:, 2*i), 'Color', 'blue');
    hold all
end
%% PLOT THE TERMINAL VELOCITY VECTORS
hold all
for i = 1:N
    h = quiver(sol(end, 2*i-1), sol(end, 2*i), 5*vT(i, 1), 5*vT(i, 2),'filled');
    h.MaxHeadSize = 1;
    set(h,'linewidth',1);
    set(h,'color',[1,0,0]);
end
% %% PLOT THE INITIAL VELOCITY VECTORS
% hold all
% for i = 1:N
%     h = quiver(initial_x0(i, 1), initial_x0(i, 2), 1*initial_v0(i, 1), 1*initial_v0(i, 2),'filled');
%     h.MaxHeadSize = 1;
%     set(h,'linewidth',1);
%     set(h,'color',[1,0,0]);
% end
% title('evolution');



% %% PLOT THE CONTROLS
% % d = 1
% figure
% for i = 1:N
%     plot(t(1:end-ndT), control(2*i-1, :, 1));
%     hold all
% end
% title('controls d = 1');
% 
% % d = 2
% figure
% for i = 1:N
%     plot(t(1:end-ndT), control(2*i, :, 1));
%     hold all
% end
% title('controls d = 2');

% 
% %% PLOT X
% figure
% for k = 1:length(t)
%     x = reshape(sol(k, 1 : N*d), [d, N])';
%     YX(k) =  B(x, x, N);
% end
% plot(t, YX);
% title('X(t)');
% 
% 
% 
% %% PLOT E
% figure
% for k = 1:length(t)
%     x = reshape(sol(k, 1 : N*d), [d, N])';
%     v = reshape(sol(k, N*d+1 : 2*N*d), [d, N])';
%     YE(k) =  E(x, v, N);
% end
% plot(t, YE);
% title('E(t)');

%% PLOT THE GRAPH
figure
% hold all
[X, Y] = gplot(Adjc, xT);
plot(X, Y, '-o', 'Color', 'b');
title('connectivity graph');

hold all
Adjc0 = get_adjacency(initial_x0, N, R, Rh, zeros(N));
[X, Y] = gplot(Adjc, initial_x0);
plot(X, Y, '-o', 'Color', 'b');

%% PLOT THE TERMINAL VELOCITY VECTORS
hold all
for i = 1:N
    h = quiver(sol(end, 2*i-1), sol(end, 2*i), 5*vT(i, 1), 5*vT(i, 2),'filled');
    h.MaxHeadSize = 1;
    set(h,'linewidth',1);
    set(h,'color',[1,0,0]);
end
title('evolution');
%% PLOT THE INITIAL VELOCITY VECTORS
hold all
for i = 1:N
    h = quiver(initial_x0(i, 1), initial_x0(i, 2), 1*initial_v0(i, 1), 1*initial_v0(i, 2),'filled');
    h.MaxHeadSize = 1;
    set(h,'linewidth',1);
    set(h,'color',[0,0,0]);
end
