% %% GET ENDTIME VALUES
% [xT, vT, zT, uT] = convert(solx(:, end), solu(:, end, 1), N, d);
% 
% 
% %% NORM of the SYSTEM VELOCITY at the end-time
% normv = norm(solx(N*d+1:2*N*d, end))
% 
% 
% 
% %% PLOT THE LYAPUNOV FUNCTION
% figure
% for k = 1:length(t)
% %     x = reshape(sol(k, 1 : N*d), [d, N])';
%     v = reshape(sol(k, N*d+1 : 2*N*d), [d, N])';
%     YV(k) =  B(v, N, d);
% end
% plot(t, YV);
% title('V(t) = 1/2N^2  sumij||vi -vj ||^2');
% 
% 
% 
% 
% %% PLOT TRAJECTORIES
% figure
% for i = 1:N
%     plot(sol(:, 2*i-1), sol(:, 2*i));
%     hold all
% end
% title('evolution');
% 
% 
% %% PLOT THE CONTROLS
% % d = 1
% figure
% for i = 1:N
%     plot(t(1:end-1), solu(2*i-1, :, 1));
%     hold all
% end
% title('controls d = 1');
% 
% % d = 2
% figure
% for i = 1:N
%     plot(t(1:end-1), solu(2*i, :, 1));
%     hold all
% end
% title('controls d = 2');
% 
% 
% %% PLOT X
% figure
% for k = 1:length(t)
%     x = reshape(sol(k, 1 : N*d), [d, N])';
%     YX(k) =  B(x, N, d);
% end
% plot(t, YX);
% title('X(t)');
% 
% 

display_results(solution, control, N, d, T, x0, v0);
