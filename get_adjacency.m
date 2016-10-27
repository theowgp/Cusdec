function res = get_adjacency(x, N, R, Rh, Adjc)
res = Adjc;

for i = 1:N
    for j = 1:N
        if i ~= j
            if norm(x(i, :) - x(j, :)) < Rh 
                res(i, j) = 1;
            else
                if norm(x(i, :) - x(j, :)) >= R
                    res(i, j) = 0;
                end
            end
        end
    end
end

end

