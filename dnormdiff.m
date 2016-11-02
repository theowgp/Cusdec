function res = dnormdiff(w, i, j, k, d)
            res = zeros(1, d);
            if i ~= j
                if k == i
                    res = (w(k, :) - w(j, :))/norm(w(k, :) - w(j, :));
                else
                    if k == j
                        res = -(w(i, :) - w(k, :))/norm(w(i, :) - w(k, :));
                    end
                end
            end
end