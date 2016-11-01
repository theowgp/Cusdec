function res = dnormdiff(w, i, j, k, d)
            res = zeros(1, d);
            if i ~= j
                if k == i
                    res = w(k, :)/norm(w(k, :) - w(j, :));
                else
                    if k == j
                        res = -w(k, :)/norm(w(i, :) - w(k, :));
                    end
                end
            end
end