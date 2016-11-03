function res = ComputeBetaK(g, gnext, drct, mesh)

y = gnext - g;


temp = spsolu(y, drct, mesh);

sigma = y - 2*drct*normsolu(y, mesh)^2/temp;

res = spsolu(sigma, gnext, mesh)/temp;

end