function res = mean(v, N, d)
res = zeros(1, d);

for i = 1:N
    res = res+  v(i, :);
end

res = res/N;
end

