%FUNCTION HANDLES - Rosenbrock function
f_Rosenbrock = @(x) 100*(x(2) - x(1)^2)^2 + (1-x(1))^2;
grad_Rosenbrock = @(x) [400*x(1)^3 - 400*x(1)*x(2) + 2*x(1) - 2; 200*(x(2) - x(1)^2)];
Hess_Rosenbrock = @(x) [1200*x(1)^2 - 400*x(2)+2,     -400*x(1); -400*x(1),     200];

%Starting points
x0_1 = [1.2; 1.2];
x0_2 = [-1.2; 1];

%Parameters
kmax = 1000;
rho = 0.5;                          
c = 10^(-4);
h = sqrt(eps);                     
btmax = 100;


% For pcg (if used)
tolgrad = 10^-12;
pcg_maxit = 100 ;


%% STANDARD NEWTON METHOD:--------------------------------------------------------------------------------------------------
%{
%EXACT - EXACT
FDgrad = '';
FDHess = '';

tStart_1 = cputime;
[xk_1, fk_1, gradfk_norm_1, k_1, non_positive_1, xseq_1, ~, ~, ~, exact_time_1] = newton_backtrack(x0_1, f_Rosenbrock, grad_Rosenbrock, ...
                                                  Hess_Rosenbrock, kmax, tolgrad, c, rho, btmax, FDgrad, FDHess, h, pcg_maxit);
tEnd_1 = cputime - tStart_1
                                              
tStart_2 = cputime;
[xk_2, fk_2, gradfk_norm_2, k_2, non_positive_2, xseq_2,  ~, ~, ~, exact_time_2] = newton_backtrack(x0_2, f_Rosenbrock, grad_Rosenbrock, ...
                                                  Hess_Rosenbrock, kmax, tolgrad, c, rho, btmax, FDgrad, FDHess, h, pcg_maxit);
tEnd_2 = cputime - tStart_2


%}
%% STANDARD NEWTON METHOD: finite difference derivatives approximation------------------------------------------------------------------------

%{                            
%EXACT and JACOBIAN CENTERED
FDgrad = '';
FDHessf = 'Jc';

tStart_1_Jc = cputime;
[xk_1_Jc, fk_1_Jc, gradfk_norm_1_Jc, k_1_Jc, non_positive_1_Jc, xseq_1_Jc, ~, ~, ~, exact_time_1_Jc] = newton_backtrack(x0_1, f_Rosenbrock, grad_Rosenbrock, ...
                                                  Hess_Rosenbrock, kmax, tolgrad, c, rho, btmax, FDgrad, FDHess, h, pcg_maxit)
tEnd_1_Jc = cputime - tStart_1_Jc

tStart_2_Jc = cputime;
[xk_2_Jc, fk_2_Jc, gradfk_norm_2_Jc, k_2_Jc, non_positive_2_Jc, xseq_2_Jc, ~, ~, ~, exact_time_2_Jc] = newton_backtrack(x0_2, f_Rosenbrock, grad_Rosenbrock, ...
                                                  Hess_Rosenbrock, kmax, tolgrad, c, rho, btmax, FDgrad, FDHess, h, pcg_maxit) 
tEnd_2_Jc = cputime - tStart_2_Jc
%}

%{
%CENTERED - FORWARD
FDgrad = 'c';
FDHessf = 'fw';

tStart_1_cf = cputime;
[xk_1_cf, fk_1_cf, gradfk_norm_1_cf, k_1_cf, non_positive_1_cf, xseq_1_cf, ~, ~, ~, exact_time_1_cf] = newton_backtrack(x0_1, f_Rosenbrock, grad_Rosenbrock, ...
                                                  Hess_Rosenbrock, kmax, tolgrad, c, rho, btmax, FDgrad, FDHess, h, pcg_maxit)
tEnd_1_cf = cputime - tStart_1_cf
                                              
tStart_2_cf = cputime;
[xk_2_cf, fk_2_cf, gradfk_norm_2_cf, k_2_cf, non_positive_2_cf, xseq_2_cf, ~, ~, ~, exact_time_2_cf] = newton_backtrack(x0_2, f_Rosenbrock, grad_Rosenbrock, ...
                                                  Hess_Rosenbrock, kmax, tolgrad, c, rho, btmax, FDgrad, FDHess, h, pcg_maxit)
tEnd_2_cf = cputime - tStart_1_cf                                    
%}      

%MATRIX FREE IMPLEMENTATION - We use pcg            
%{
FDgrad = '';
FDHess = 'MF';

tStart_1_MF = cputime;
[xk_1_MF, fk_1_MF, gradfk_norm_1_MF, k_1_MF, non_positive_1_MF, xseq_1_MF] = newton_backtrack(x0_1, f_Rosenbrock, grad_Rosenbrock, ...
                                                   Hess_Rosenbrock, kmax, tolgrad, c, rho, btmax, FDgrad, FDHess, h, pcg_maxit)
tEnd_1_MF = cputime - tStart_1_MF

tStart_2_MF = cputime;
[xk_2_MF, fk_2_MF, gradfk_norm_2_MF, k_2_MF, non_positive_2_MF, xseq_2_MF] = newton_backtrack(x0_2, f_Rosenbrock, grad_Rosenbrock, ...
                                                   Hess_Rosenbrock, kmax, tolgrad, c, rho, btmax, FDgrad, FDHess, h, pcg_maxit)
tEnd_2_MF = cputime - tStart_2_MF
%}


%% EXACT RESULTS
%{
xmin = [1; 1];
fmin = f_Rosenbrock(xmin);

fprintf('Exact value of xmin:  (%2.5f, %2.5f) \n', xmin(1), xmin(2));
fprintf('Exact value of the function in xmin: %2.5f \n', fmin);
%}


%% PLOTS---------------------------------------------------------------------------------------------------------------------------------------
%{
lowbound = 0.80;
upbound = 1.50;

plot_results(f_Rosenbrock, x0_1, xseq_1, lowbound, upbound);
plot_results(f_Rosenbrock, x0_1, xseq_1_Jc, lowbound, upbound);
plot_results(f_Rosenbrock, x0_1, xseq_1_cf, lowbound, upbound);
plot_results(f_Rosenbrock, x0_1, xseq_1_MF, lowbound, upbound);
%}
%{
lowbound = -1.5;
upbound = 1.6;

plot_results(f_Rosenbrock, x0_2, xseq_2, lowbound, upbound);
plot_results(f_Rosenbrock, x0_2, xseq_2_Jc, lowbound, upbound);
plot_results(f_Rosenbrock, x0_2, xseq_2_cf, lowbound, upbound);
plot_results(f_Rosenbrock, x0_2, xseq_2_MF, lowbound, upbound);


%}


%% TEST 
%
FDgrad = '';
FDHess = '';
[xsol_1, gradsol_1, ksol_1, tsol_1] = test_solver_kmax(f_Rosenbrock, grad_Rosenbrock, Hess_Rosenbrock, x0_1, FDgrad, FDHess, 2)

[xsol_2, gradsol_2, ksol_2, tsol_2] = test_solver_kmax(f_Rosenbrock, grad_Rosenbrock, Hess_Rosenbrock, x0_2, FDgrad, FDHess, 2)
%}

%
FDgrad = '';
FDHess = 'Jc';

[xsol_1_JC, gradsol_1_JC, ksol_1_JC, tsol_1_JC] = test_solver_kmax(f_Rosenbrock, grad_Rosenbrock, Hess_Rosenbrock, x0_1, FDgrad, FDHess, 2)

[xsol_2_JC, gradsol_2_JC, ksol_2_JC, tsol_2_JC] = test_solver_kmax(f_Rosenbrock, grad_Rosenbrock, Hess_Rosenbrock, x0_2, FDgrad, FDHess, 2)
%}

%
FDgrad = 'c';
FDHess = 'fw';

[xsol_1_cf, gradsol_1_cf, ksol_1_cf, tsol_1_cf] = test_solver_kmax(f_Rosenbrock, grad_Rosenbrock, Hess_Rosenbrock, x0_1, FDgrad, FDHess, 2)

[xsol_2_cf, gradsol_2_cf, ksol_2_cf, tsol_2_cf] =  test_solver_kmax(f_Rosenbrock, grad_Rosenbrock, Hess_Rosenbrock, x0_2,FDgrad, FDHess, 2)
%}

%
FDgrad = '';
FDHess = 'MF';
[xsol_1_MF, gradsol_1_MF, ksol_1_MF, tsol_1_MF] = test_solver_kmax(f_Rosenbrock, grad_Rosenbrock, Hess_Rosenbrock, x0_1, FDgrad, FDHess, 2)
FDgrad = '';
FDHess = 'MF';
[xsol_2_MF, gradsol_2_MF, ksol_2_MF, tsol_2_MF] = test_solver_kmax(f_Rosenbrock, grad_Rosenbrock, Hess_Rosenbrock, x0_2, FDgrad, FDHess, 2)
%}
















