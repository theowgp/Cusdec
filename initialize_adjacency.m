function res = initialize_adjacency(x, N, R)
res = zeros(N);

for i = 1:N
    for j = 1:N
        if i ~= j
            if norm(x(i, :) - x(j, :)) < R - 0.1
                res(i, j) = 1;
            end
        end
    end
end

end
