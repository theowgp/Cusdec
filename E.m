function res = E(x, v, N, d)

V = B(v, N, d);
X = B(x, N, d);

res = sqrt(V);


temp = (1/sqrt(2*N)) * (pi/2 - atan(sqrt(2*N*X)));

res = res-  temp;


end

