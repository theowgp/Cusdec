function [xT, vT] =  display_results(solx, solu, N, d, T, x0, v0, Adjc, R)


[m, n] = size(solx) ;
mesh = Mesh(T, n-1);
t = mesh.t;

sol = solx';

%% GET ENDTIME VALUES
[xT, vT, zT, uT] = convert(solx(:, end), solu(:, end, 1), N, N, d);


%% NORM of the SYSTEM VELOCITY at the end-time
normv = norm(solx(N*d+1:2*N*d, end))



%% PLOT THE LYAPUNOV FUNCTION
figure
for k = 1:length(t)
    [x, v, z] = convert_state(sol(k, :), N, d);
    YV(k) =  B(v, v, N);
end
plot(t, YV);





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
% d = 2
figure
for i = 1:N
    plot(tu(1:end), solu(2*i, :, 1));
    hold all
end

%% PLOT THE NORM OF CONTROLS
figure
tempu = 0;
for k = 1:length(tu)
     u = convert_control(solu(:, k, 1), N, d);
     YU(k) = norm(u)^2;
     tempu = tempu + norm(u)^2;
end
ControlEnergy = tempu * meshu.h;
plot(tu, YU);





%% PLOT THE GRAPH
figure
% hold all
[X, Y] = gplot(Adjc, xT);
plot(X, Y, '-o', 'Color', 'b');
% title('connectivity graph');

hold all
Adjc0 = get_adjacency(x0, N, R, zeros(N));
[X, Y] = gplot(Adjc0, x0);
plot(X, Y, '-o', 'Color', 'b');

%% PLOT THE TERMINAL VELOCITY VECTORS
hold all
for i = 1:N
    h = quiver(sol(end, 2*i-1), sol(end, 2*i), 5*vT(i, 1), 5*vT(i, 2),'filled');
    h.MaxHeadSize = 1;
    set(h,'linewidth',1);
    set(h,'color',[1,0,0]);
end
% title('graph evolution');
%% PLOT THE INITIAL VELOCITY VECTORS
hold all
for i = 1:N
    h = quiver(x0(i, 1), x0(i, 2), 0.5*v0(i, 1), 0.5*v0(i, 2),'filled');
    h.MaxHeadSize = 0.5;
    set(h,'linewidth',1);
    set(h,'color',[0,0,0]);
end


BvT = B(vT, vT, N)

ET = E(vT, vT, N)

vT


ControlEnergy


g0 = graph(Adjc0);
gT = graph(Adjc);

figure
plot(g0);
figure
plot(gT);








BvT = B(vT, vT, d)

ET = E(x, v, N)


vT







