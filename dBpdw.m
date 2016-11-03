function res = dBpdw(w, p, k, N, d)

res = zeros(1, d);

if k == p
    for j = 1:N
        if j ~= k
            res = res+  2 * norm(w(k, :) - w(j, :)) * dnormdiff(w, k, j, k, d);
        end
    end
else
    res = res+  2 * norm(w(p, :) - w(k, :)) * dnormdiff(w, p, k, k, d);
end

res = res/N;
    

end

