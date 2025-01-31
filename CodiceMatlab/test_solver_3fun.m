function [] = test_solver(f, grad, Hess, x0_1, x0_2, FDgrad, FDHess)

%Fixed parameters
rho = 0.5;                          
h = sqrt(eps);                     

%Variable parameters
kmax = [100, 10000, 1000000];               %i
c = [10^(-2), 10^(-4)];                     %j
tolgrad = [10^(-4), 10^(-8), 10^(-12)];     %q
btmax = [10, 100, 1000];                    %u
pcg_maxit = [10, 100, 1000];                %v

j = 2;
q = 3;
u = 1;
v = 1;

xsol_1 = zeros(3, 2);           %per righe
xsol_2 = zeros(3,2);            %per colonne

gradsol_1 = zeros(3),
gradsol_2 = zeros(3);

ksol_1 = [0, 0, 0];
ksol_2 = [0, 0, 0];

tsol_1 = [0, 0, 0];
tsol_2 = [0, 0, 0];

for i = 1:3
    
    tStart_1= cputime;
    [xk_1, fk_1, gradfk_norm_1, k_1] = newton_backtrack(x0_1, f, grad, Hess, kmax(i),tolgrad(q), ...
                                                                c(j), rho, btmax(u), FDgrad, FDHess, h, pcg_maxit(v));
    tEnd_1 = cputime - tStart_1;

    tStart_2 = cputime;
    [xk_2, fk_2, gradfk_norm_2, k_2] = newton_backtrack(x0_2, f, grad, Hess, kmax (i), tolgrad(q), ....
                                                                 c(j), rho, btmax(u), FDgrad, FDHess, h, pcg_maxit(v));
    tEnd_2 = cputime - tStart_2;
    
    
    xsol_1(i, :) = xk_1;
    xsol_2(i, :) = xk_2;
    gradsol_1(i, :) = gradfk_norm_1;            %per righe
    gradsol_2(i, :) = gradfk_norm_2;            %per righe
    ksol_1(i) = k_1;
    ksol_2(i) = k_2;
    tsol_1(i)= tEnd_1;
    tsol_2(i) = tEnd_2;
    
    i = i+1;
end


end