function [solx, solu] = NCG(rk, objective, mesh, solu, eps, sigma, limitLS, limitA)

g =  rk.g_u(solu);
drct = - g;


Phi0 = 9999999999999999;
[solx, soly] = rk.solve_forward_equation(solu);
Phi1 = objective.phi(solx(:, end));


kLS = 0;
while  kLS<limitLS && Phi1 < Phi0
    [step, kA] = DetermineStepSize(rk, objective, mesh, solu, g, drct, sigma, limitA);
    solu = solu + step * drct;

    
    gnext = rk.g_u(solu);
    beta = ComputeBetaK(g, gnext, drct, mesh);
    drct = - gnext + beta*drct;
    g = gnext;
    
    
    
    kLS = kLS+1;
    
    
    
    Gradient = normsolu(g, mesh);
    
    Phi0 = Phi1;
    
    [solx, soly] = rk.solve_forward_equation(solu);
    Phi1 = objective.phi(solx(:, end));
    
    disp([kLS, kA, Phi1, Gradient]); 
end



end
