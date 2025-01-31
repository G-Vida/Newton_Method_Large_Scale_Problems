function [xsol, gradsol, ksol, tsol] = test_solver_kmax(f, grad, Hess, x0, FDgrad, FDHess, n)

%Fixed parameters
rho = 0.5;                          
h = sqrt(eps);                     
c = 10^(-4);                     
tolgrad = 10^(-12);     
btmax = 100;                    
pcg_maxit = 100;                


%Variable parameters
kmax = [10000, 100000, 1000000];               %i


xsol = zeros(3, n);  
gradsol = zeros(3, 1),
ksol = [0, 0, 0];
tsol_1 = [0, 0, 0];


% Varying kmax
for i = 1:3
    
    tStart= cputime;
    [xk, ~, gradfk_norm, k] = newton_backtrack(x0, f, grad, Hess, kmax(i),tolgrad, ...
                                                                c, rho, btmax, FDgrad, FDHess, h, pcg_maxit);
    tEnd = cputime - tStart;

    
   
    xsol(i, :) = xk;
    gradsol(i) = gradfk_norm;
    ksol(i) = k;
    tsol(i)= tEnd;
    
    i = i+1;
end


end