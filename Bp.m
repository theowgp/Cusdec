function res = Bp(w, p, N)
res = 0;

for j = 1:N
    if j ~= p
        res = res+  norm(w(p, :) - w(j, :))^2;
    end
end

res = res/N;
end

