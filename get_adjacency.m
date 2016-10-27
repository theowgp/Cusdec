function res = get_adjacency(x, N, R)
res = zeros(N);

for i = 1:N
    for j = 1:N
        if norm(x(i, :) - x(j, :)) <= R && i ~= j
            res(i, j) = 1;
        end
    end
end

end

