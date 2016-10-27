function w0 = set_w0(arg0, A, N, i)

w0(1, :) = arg0(i, :);

s = 2;
for j = 1:N
    if A(i, j) == 1 && i ~= j
        w0(s, :) = arg0(j, :);
        s =  s + 1;
    end
end

end

