function v = initv( N, d, maxv )
v = zeros(N, d);

% a=2;
% b=-1;
% 
% for i=1:d
%     t = rand(N, 1);
%     v(:, i) = (a*t + b) * maxv;
% end


% v =  0.5*ones(N, d);

v = [
    -1 0;
     0 1;
     1 1
    ];
end
