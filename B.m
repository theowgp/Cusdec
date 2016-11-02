function res = B(v, N, d)

vbar = mean(v, N, d);

temp = 0;
for i = 1:N
    temp = temp+  norm(v(i, :) - vbar)^2;
%     temp = temp+  norm(v(i, :) - vbar);
end

res = temp/N;
end
