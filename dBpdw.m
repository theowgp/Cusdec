function res = dBpdw(w, p, k, N, d)

res = zeros(1, d);

for j = 1:N
    if j ~= p
        res = res+  2 * norm(w(p, :) - w(j, :)) * dnormdiff(w, p, j, k, d);
    end
end

res = res/N;
    

end

