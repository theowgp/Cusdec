function res = dnormvmeandiff(v, i, k, N, d)

vbar = mean(v, N, d);

if k == i
    temp = v(k, :) - vbar;
    res = temp * (1 - 1/N) / norm(temp);
else
    temp = v(i, :) - vbar;
    res = temp * (- 1/N) / norm(temp);
end

end

