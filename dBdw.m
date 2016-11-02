function res = dBdw(v, k, N, d)

vbar = mean(v, N, d);

temp = zeros(1, d);
for i = 1:N
    temp = temp+    2 * norm(v(i, :) - vbar) * dnormvmeandiff(v, i, k, N, d);
%     temp = temp+    dnormvmeandiff(v, i, k, N, d);
end

res = temp/N;
end
