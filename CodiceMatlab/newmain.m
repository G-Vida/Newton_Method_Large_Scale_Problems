%FUNCTION HANDLES - Rosenbrock function
f_Rosenbrock = @(x) 100*(x(2) - x(1)^2)^2 + (1-x(1))^2;
grad_Rosenbrock = @(x) [400*x(1)^3 - 400*x(1)*x(2) + 2*x(1) - 2; 200*(x(2) - x(1)^2)];
Hess_Rosenbrock = @(x) [1200*x(1)^2 - 400*x(2)+2,     -400*x(1); -400*x(1),     200];

%Starting point
x0_1 = [1.2; 1.2];
x0_2 = [-1.2; 1];


%Parameters
kmax = 100;
rho = 0.5;                          
c = 10^(-4);
h = sqrt(eps);                     
tolgrad = 10^-12;
btmax = 50;
FDgrad = 'fw';
FDHess = 'fw';
pcg_maxit = 50;





%% STANDARD NEWTON METHOD: exact derivatives--------------------------------------------------------------------------------------------------
FDgrad = '';
FDHess = '';

tStart_1 = cputime;
[xk_1, fk_1, gradfk_norm_1, k_1, xseq_1, btseq_1, non_positive_1] = newton_backtrack(x0_1, f_Rosenbrock, grad_Rosenbrock, ...
                                                  Hess_Rosenbrock, kmax, tolgrad, c, rho, btmax, FDgrad, FDHess, h, pcg_maxit);
tEnd_1 = cputime - tStart_1;
                                              
tStart_2 = cputime;
[xk_2, fk_2, gradfk_norm_2, k_2, xseq_2, btseq_2, non_positive_2] = newton_backtrack(x0_2, f_Rosenbrock, grad_Rosenbrock, ...
                                                  Hess_Rosenbrock, kmax, tolgrad, c, rho, btmax, FDgrad, FDHess, h, pcg_maxit);
tEnd_2 = cputime - tStart_2;




%% STANDARD NEWTON METHOD: finite difference derivatives approximation------------------------------------------------------------------------

%FORWARD and FORWARD
FDgrad = 'fw';
FDHessf = 'fw';

[xk_1_ff, fk_1_ff, gradfk_norm_1_ff, k_1_ff, xseq_1_ff, btseq_1_ff, non_positive_1_ff] = newton_backtrack(x0_1, f_Rosenbrock, grad_Rosenbrock, ...
                                                  Hess_Rosenbrock, kmax, tolgrad, c, rho, btmax, FDgrad, FDHess, h, pcg_maxit);
[xk_2_ff, fk_2_ff, gradfk_norm_2_ff, k_2_ff, xseq_2_ff, btseq_2_ff, non_positive_2_ff] = newton_backtrack(x0_2, f_Rosenbrock, grad_Rosenbrock, ...
                                                  Hess_Rosenbrock, kmax, tolgrad, c, rho, btmax, FDgrad, FDHess, h, pcg_maxit); 

                                              
%CENTERED and CENTERED
FDgrad = 'c';
FDHessf = 'c';

tStart_1_cc = cputime;
[xk_1_cc, fk_1_cc, gradfk_norm_1_cc, k_1_cc, xseq_1_cc, btseq_1_cc, non_positive_1_cc] = newton_backtrack(x0_1, f_Rosenbrock, grad_Rosenbrock, ...
                                                  Hess_Rosenbrock, kmax, tolgrad, c, rho, btmax, FDgrad, FDHess, h, pcg_maxit);
tEnd_1_cc = cputime - tStart_1_cc;

tStart_2_cc = cputime;
[xk_2_cc, fk_2_cc, gradfk_norm_2_cc, k_2_cc, xseq_2_cc, btseq_2_cc, non_positive_2_cc] = newton_backtrack(x0_2, f_Rosenbrock, grad_Rosenbrock, ...
                                                  Hess_Rosenbrock, kmax, tolgrad, c, rho, btmax, FDgrad, FDHess, h, pcg_maxit); 
tEnd_2_cc = cputime - tStart_2_cc;



%FORWARD and FORWARD-Jacobian of the Gradient
FDgrad = 'fw';
FDHessf = 'Jfw';

[xk_1_fJf, fk_1_fJf, gradfk_norm_1_fJf, k_1_fJf, xseq_1_fJf, btseq_1_fJf, non_positive_1_fJf] = newton_backtrack(x0_1, f_Rosenbrock, grad_Rosenbrock, ...
                                                  Hess_Rosenbrock, kmax, tolgrad, c, rho, btmax, FDgrad, FDHess, h, pcg_maxit);
[xk_2_fJf, fk_2_fJf, gradfk_norm_2_fJf, k_2_fJf, xseq_2_fJf, btseq_2_fJf, non_positive_2_fJf] = newton_backtrack(x0_2, f_Rosenbrock, grad_Rosenbrock, ...
                                                  Hess_Rosenbrock, kmax, tolgrad, c, rho, btmax, FDgrad, FDHess, h, pcg_maxit); 

                                              
                                              
%CENTERED and CENTERED-Jacobian of the Gradient
FDgrad = 'c';
FDHessf = 'Jc';

[xk_1_cJc, fk_1_cJc, gradfk_norm_1_cJc, k_1_cJc, xseq_1_cJc, btseq_1_cJc, non_positive_1_cJc] = newton_backtrack(x0_1, f_Rosenbrock, grad_Rosenbrock, ...
                                                    Hess_Rosenbrock, kmax, tolgrad, c, rho, btmax, FDgrad, FDHess, h, pcg_maxit);
[xk_2_cJc, fk_2_cJc, gradfk_norm_2_cJc, k_2_cJc, xseq_2_cJc, btseq_2_cJc, non_positive_2_cJc] = newton_backtrack(x0_2, f_Rosenbrock, grad_Rosenbrock, ...
                                                    Hess_Rosenbrock, kmax, tolgrad, c, rho, btmax, FDgrad, FDHess, h, pcg_maxit); 

                                                
                                              
%MATRIX FREE IMPLEMENTATION             
FDgrad = '';
FDHess = 'MF';

tStart_1_MF = cputime;
[xk_1_MF, fk_1_MF, gradfk_norm_1_MF, k_1_MF, xseq_1_MF, btseq_1_MF, non_positive_1_MF] = newton_backtrack(x0_1, f_Rosenbrock, grad_Rosenbrock, ...
                                                   Hess_Rosenbrock, kmax, tolgrad, c, rho, btmax, FDgrad, FDHess, h, pcg_maxit);
tEnd_1_MF = cputime - tStart_1_MF;

tStart_2_MF = cputime;
[xk_2_MF, fk_2_MF, gradfk_norm_2_MF, k_2_MF, xseq_2_MF, btseq_2_MF, non_positive_2_MF] = newton_backtrack(x0_2, f_Rosenbrock, grad_Rosenbrock, ...
                                                   Hess_Rosenbrock, kmax, tolgrad, c, rho, btmax, FDgrad, FDHess, h, pcg_maxit);
tEnd_2_MF = cputime - tStart_2_MF;





%% PLOTS---------------------------------------------------------------------------------------------------------------------------------------
lowbound = 0.8;
upbound = 1.5;

plot_results(f_Rosenbrock, x0_1, xseq_1, lowbound, upbound);
%plot_results(f_Rosenbrock, x0_1, xseq_1_ff, lowbound, upbound);
plot_results(f_Rosenbrock, x0_1, xseq_1_cc, lowbound, upbound);
%plot_results(f_Rosenbrock, x0_1, xseq_1_fJf, lowbound, upbound);
%plot_results(f_Rosenbrock, x0_1, xseq_1_cJc, lowbound, upbound);
plot_results(f_Rosenbrock, x0_1, xseq_1_MF, lowbound, upbound);




lowbound = -1.4;
upbound = 1.9;

plot_results(f_Rosenbrock, x0_2, xseq_2, lowbound, upbound);
%plot_results(f_Rosenbrock, x0_2, xseq_2_ff, lowbound, upbound);
plot_results(f_Rosenbrock, x0_2, xseq_2_cc, lowbound, upbound);
%plot_results(f_Rosenbrock, x0_2, xseq_2_fJf, lowbound, upbound);
%plot_results(f_Rosenbrock, x0_2, xseq_2_cJc, lowbound, upbound);
plot_results(f_Rosenbrock, x0_2, xseq_2_MF, lowbound, upbound);


%% RESULTS - TIME OF CONVERGENCE
disp(['Result of exact method starting from point 1:  ', num2str(xk_1)])


%% RESULTS - TIME OF CONVERGENCE 
%I tempi restano significativi anche se conprendenti il backtracking, 
%perchè l'operazione più costosa computaizonalmente è il pcg (risoluzione sistmea lineare)

disp(['CPU time (in s) of exact method, point 1:  ', num2str(tEnd_1)])
disp(['CPU time (in s) of exact method, point 2:  ', num2str(tEnd_2)])

disp(['CPU time (in s) of centered-centered method, point 1:  ', num2str(tEnd_2_cc)])
disp(['CPU time (in s) of centered-centered method, point 2:  ', num2str(tEnd_2_cc)])

disp(['CPU time (in s) of matrix free method, point 1:  ', num2str(tEnd_1_MF)])
disp(['CPU time (in s) of matrix free method, point 2:  ', num2str(tEnd_2_MF)])





xmin = [1; 1];

%fmin = f_Rosenbrock(xmin)


%f0_1 = f_Rosenbrock(x0_1)
%f0_2 = f_Rosenbrock(x0_2)
%xk_1
%xk_2
%fk_1 = f_Rosenbrock(xk_1)
%fk_2 = f_Rosenbrock(xk_2)

%xmin
%fmin = f_Rosenbrock(xmin)








