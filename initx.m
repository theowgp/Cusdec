function x = initx(N, d, maxx, R )
x = zeros(N, d);

% % random
% for i=1:d
%     t = rand(N, 1);
%     x(:, i) = t * maxx;
% end


% % uniform
% for i = 2:N
%     x(i, 1) = x(i-1, 1) + 1;
% end

% x = [
%     -2 0
%     0 0
%     2 0
%     ];

% x = [
%      1  1  
%     -1  0
%      0  1
%      1  0
%      0 -1
%     ];


% x = [
%      0  0  
%      1  1
%      2  0
%      3  1
%      4  0
%      5  1
%      6  0
%      7  1
%      8  0
%      9  1
%      
%     ];



x = [
    -2 0
    0 2
    2 0
    ];


% 
% % uniform circle
% eps = 0.3;
% 
% radius = N*R/(2*pi) - eps;
% t = linspace(0, 2*pi, N+1);
% 
% for i=1:N
%     x(i, 1) = radius *cos(t(i));
%     x(i, 2) = radius *sin(t(i));
% end


end