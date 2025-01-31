function [fsol, gradsol, ksol, tsol] = test_solver_random_points(f, grad, Hess, FDgrad, FDHess, kmax, n)


%Fixed parameters
rho = 0.5;                          
h = sqrt(eps);                     
c = 10^(-4);                     
tolgrad = 10^(-12);     
btmax = 100;                    
pcg_maxit = 100; 

%Variable parameters


%Random starts
rng(11);
xr0_1 = randi([-4, 5], 1, n)*0.75;
xr0_1 = xr0_1';

rng(15);
xr0_2 = randi([-4, 5], 1, n)*0.75;
xr0_2 = xr0_2';

rng(23);
xr0_3 = randi([-4, 5], 1, n)*0.75;
xr0_3 = xr0_3';

xr = [xr0_1, xr0_2, xr0_3];
 


fsol = [0, 0, 0];
gradsol = zeros(3, 1);
ksol = [0, 0, 0];
tsol = [0, 0, 0];
  

for i = 1:3
        
    tStart= cputime;
    [~, fk, gradfk_norm, k, ~, ~, ~, ~] = newton_backtrack(xr(:, i), f, grad, Hess, kmax,tolgrad, ...
                                                                c, rho, btmax, FDgrad, FDHess, h, pcg_maxit);
    tEnd = cputime - tStart;

    
   
    fsol (i) = fk;
    gradsol(i) = gradfk_norm;
    ksol(i) = k;
    tsol(i)= tEnd;
    
    i = i+1;

end

end