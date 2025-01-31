%% Test for high dimensions, starting from suggested points-----------------------------------------------------------------------------

%
%Parameters
kmax = 10000;
rho = 0.5;
c = 10^(-4);
h = sqrt(eps);
btmax = 100;

% For pcg (when used)
tolgrad = 10^-12;
pcg_maxit = 100;

n=5;


%% PROBLEM 1 - Chained Rosenbrock function                             
%{

FDgrad = 'c';
FDHess = 'fw';

f_Chained_Rosenbrock = @fun_f_Chained_Rosenbrock;
f_gradient_Chained_Rosenbrock = @grad_f_Chained_Rosenbrock;
f_Hessian_Chained_Rosenbrock = @Hess_f_Chained_Rosenbrock;


%suggested start
x0 = [-1.2; 1.0; -1.2; 1.0; -1.2];

[xk_c, fk_c, gradfk_norm_c, k_c, non_positive_c] = newton_backtrack(x0, f_Chained_Rosenbrock, f_gradient_Chained_Rosenbrock , ...
                                    f_Hessian_Chained_Rosenbrock, kmax, tolgrad, c, rho, btmax, FDgrad, FDHess, h, pcg_maxit)                        
%}

   
%% PROBLEM 16 - Banded trigonometric problem

%{

FDgrad = '';
FDHess = '';

f_Banded = @fun_f_Banded;
f_gradient_Banded = @grad_f_Banded;
f_Hessian_Banded = @Hess_f_Banded;

%suggested start
x0 = [1; 1; 1; 1; 1];

[xk_ba, fk_ba, gradfk_norm_ba, k_ba, non_positive_ba] = newton_backtrack(x0, f_Banded, f_gradient_Banded , ...
                                    f_Hessian_Banded, kmax, tolgrad, c, rho, btmax, FDgrad, FDHess, h, pcg_maxit)
%}

%% PROBLEM 31 - Brodyen tridiagonal function

%{
FDgrad = '';
FDHess = 'MF';

f_Broyden = @fun_f_Broyden;
f_gradient_Broyden = @grad_f_Broyden;
f_Hessian_Broyden = @Hess_f_Broyden;

%suggested start
x0 = [-1; -1; -1; -1; -1];

[xk_br, fk_br, gradfk_norm_br, k_br, non_positive_br] = newton_backtrack(x0, f_Broyden, f_gradient_Broyden,...
                            f_Hessian_Broyden, kmax, tolgrad, c, rho, btmax, FDgrad, FDHess, h, pcg_maxit)  
                        
%}

%}



%% Test for high dimensions, starting from suggested points

%{
FDgrad = '';
FDHess = '';

%
f_Chained_Rosenbrock = @fun_f_Chained_Rosenbrock;
f_gradient_Chained_Rosenbrock = @grad_f_Chained_Rosenbrock;
f_Hessian_Chained_Rosenbrock = @Hess_f_Chained_Rosenbrock;

[fsol_c, gradsol_c, ksol_c, tsol_c] = test_solver_n(f_Chained_Rosenbrock, f_gradient_Chained_Rosenbrock, f_Hessian_Chained_Rosenbrock, 1, FDgrad, FDHess)
%}
%{
f_Banded = @fun_f_Banded;
f_gradient_Banded = @grad_f_Banded;
f_Hessian_Banded = @Hess_f_Banded;

[fsol_ba, gradsol_ba, ksol_ba, tsol_ba] = test_solver_n( ...
                            f_Banded, f_gradient_Banded, f_Hessian_Banded, 2, FDgrad, FDHess)
%}

%{
f_Broyden = @fun_f_Broyden;
f_gradient_Broyden = @grad_f_Broyden;
f_Hessian_Broyden = @Hess_f_Broyden;

[fsol_br, gradsol_br, ksol_br, tsol_br] = test_solver_n( ...
                            f_Broyden, f_gradient_Broyden, f_Hessian_Broyden, 3, FDgrad, FDHess)
%}



%}


%% Different starting points

%{
FDgrad = '';
FDHess = 'MF';
kmax = 1000000;
n = 1000;

%{
f_Chained_Rosenbrock = @fun_f_Chained_Rosenbrock;
f_gradient_Chained_Rosenbrock = @grad_f_Chained_Rosenbrock;
f_Hessian_Chained_Rosenbrock = @Hess_f_Chained_Rosenbrock;

[fsol_c, gradsol_c, ksol_c, tsol_c] = test_solver_random_points( f_Chained_Rosenbrock,...
                             f_gradient_Chained_Rosenbrock, f_Hessian_Chained_Rosenbrock, FDgrad, FDHess, kmax, n)
%}

%{
f_Banded = @fun_f_Banded;
f_gradient_Banded = @grad_f_Banded;
f_Hessian_Banded = @Hess_f_Banded;

[fsol_ba, gradsol_ba, ksol_ba, tsol_ba] = test_solver_random_points( ...
                            f_Banded, f_gradient_Banded, f_Hessian_Banded, FDgrad, FDHess, kmax, n)
%}

%
f_Broyden = @fun_f_Broyden;
f_gradient_Broyden = @grad_f_Broyden;
f_Hessian_Broyden = @Hess_f_Broyden;

[fsol_br, gradsol_br, ksol_br, tsol_br] = test_solver_random_points( ...
                            f_Broyden, f_gradient_Broyden, f_Hessian_Broyden, FDgrad, FDHess, kmax, n)
%}



%}






%% Differences between preconditioned and not preconditioned pcg

%{
rho = 0.5;                          
h = sqrt(eps);                     
c = 10^(-4);                     
tolgrad = 10^(-12);     
btmax = 100;                    
pcg_maxit = 100; 
kmax = 500;
n = 10000;

FDgrad = '';
FDHess = '';

% Chained Rosenbrock
f_Chained_Rosenbrock = @fun_f_Chained_Rosenbrock;
f_gradient_Chained_Rosenbrock = @grad_f_Chained_Rosenbrock;
f_Hessian_Chained_Rosenbrock = @Hess_f_Chained_Rosenbrock;

   
x0 = zeros(n, 1);

[~, ~, ~, k_cr, ~, rel_diff_cr, tot_time_precond_cr, tot_time_nonprecond_cr] = newton_backtrack(x0, f_Chained_Rosenbrock, f_gradient_Chained_Rosenbrock , f_Hessian_Chained_Rosenbrock, .... 
                                 kmax, tolgrad, c, rho, btmax, FDgrad, FDHess, h, pcg_maxit)                        


figure()                     
axis_x = [0:1 :k_cr];
axis_y = [0, rel_diff_cr(1:k_cr)];
plot(axis_x, axis_y, 'blue')
title('Chained Rosenbrock function')
xlabel('Iterations')
ylabel('Relative difference')

figure()
axis_x = [0:1:k_cr];
axis_y1 = [0, tot_time_precond_cr(1:k_cr)];
axis_y2 = [0, tot_time_nonprecond_cr(1:k_cr)];
plot(axis_x, axis_y1, '--', axis_x, axis_y2, ':')
legend("Preconditioned", "Non-preconditioned")
title('Chained Rosenbrock function')
xlabel('Iterations')
ylabel('Computational time (s)')


%Banded
f_Banded = @fun_f_Banded;
f_gradient_Banded = @grad_f_Banded;
f_Hessian_Banded = @Hess_f_Banded;


x0 = ones(n, 1);

[~, ~, ~, k_ba, ~, rel_diff_ba, tot_time_precond_ba, tot_time_nonprecond_ba] = newton_backtrack(x0, f_Banded, f_gradient_Banded , f_Hessian_Banded, .... 
                                kmax, tolgrad, c, rho, btmax, FDgrad, FDHess, h, pcg_maxit);

figure()                     
axis_x = [0:1 :k_ba];
axis_y = [0, rel_diff_ba(1:k_ba)];
plot(axis_x, axis_y, 'blue')
title('Banded trigonometric problem')
xlabel('Iterations')
ylabel('Relative difference')

figure()
axis_x = [0:1:k_ba];
axis_y1 = [0, tot_time_precond_ba(1:k_ba)];
axis_y2 = [0, tot_time_nonprecond_ba(1:k_ba)];
plot(axis_x, axis_y1, '--', axis_x, axis_y2, ':')
legend("Preconditioned", "Non-preconditioned")
title('Banded trigonometric problem')
xlabel('Iterations')
ylabel('Computational time (s)')


%Broyden
f_Broyden = @fun_f_Broyden;
f_gradient_Broyden = @grad_f_Broyden;
f_Hessian_Broyden = @Hess_f_Broyden;

x0 = -ones(n, 1);

[~, ~, ~, k_br, ~, rel_diff_br, tot_time_precond_br, tot_time_nonprecond_br] = newton_backtrack(x0, f_Broyden, f_gradient_Broyden , f_Hessian_Broyden, .... 
                                 kmax, tolgrad, c, rho, btmax, FDgrad, FDHess, h, pcg_maxit);

figure()                                
axis_x = [0:1:k_br];
axis_y = [0, rel_diff_br(1:k_br)];
plot(axis_x, axis_y, 'blue')
title('Brodyen tridiagonal function')
xlabel('Iterations')
ylabel('Relative difference')

figure()
axis_x = [0:1:k_br];
axis_y1 = [0, tot_time_precond_br(1:k_br)];
axis_y2 = [0, tot_time_nonprecond_br(1:k_br)];
plot(axis_x, axis_y1, '--', axis_x, axis_y2, ':')
legend("Preconditioned", "Non-preconditioned")
title('Brodyen tridiagonal function')
xlabel('Iterations')
ylabel('Computational time (s)')

%}




%% COMPARISON WITH NELDER-MEAD
%
FDgrad = '';
FDHess = '';

kmax = 1000000;
n = [5:5:30]

%{
Mfsol_ba = zeros(length(n), 3)
Mtsol_ba =  zeros(length(n), 3)
    
for i=1: length(n)
    f_Chained_Rosenbrock = @fun_f_Chained_Rosenbrock;
    f_gradient_Chained_Rosenbrock = @grad_f_Chained_Rosenbrock;
    f_Hessian_Chained_Rosenbrock = @Hess_f_Chained_Rosenbrock;

    [fsol_c, ~, ksol_c, tsol_c] = test_solver_random_points(f_Chained_Rosenbrock,  ...
                            f_gradient_Chained_Rosenbrock, f_Hessian_Chained_Rosenbrock, FDgrad, FDHess, kmax, n(i))
    Mfsol_c(i,:) = fsol_c;
    Mtsol_c(i, :) =  tsol_c;

    f_name = "fun_f_Chained_Rosenbrock";
    [fn_c, kn_c, tn_c] = test_NelderMead_random_points(f_name, n(i), kmax)
    Mfn_c(i,:) = fn_c;
    Mtn_c(i, :) =  tn_c;
end

fig2 = figure()
ymax = max(max(Mtn_c))+ max(max(Mtn_c))/100
ylim([0, ymax])

for j=1:3
   figure()
   hold on
   
   x = [0, n];
   vet = Mtsol_c(:,j);
   y = [0, vet'];
   stairs(x, y)
   x = [0, n];
   vet = Mtn_c(:,j);
   y = [0, vet'];
   stairs(x, y)

   legend("Exact Newton Method", "Nelder Mead Method")
   title('Chained Rosenbrock function')
   xlabel('Dimension')
   ylabel('Computational time (s)')
   hold off
end



%}

%{

Mfsol_ba = zeros(length(n), 3)
Mtsol_ba =  zeros(length(n), 3)

for i=1: length(n)
    f_Banded = @fun_f_Banded;
    f_gradient_Banded = @grad_f_Banded;
    f_Hessian_Banded = @Hess_f_Banded;

    [fsol_ba, ~, ~, tsol_ba] = test_solver_random_points( ...
                            f_Banded, f_gradient_Banded, f_Hessian_Banded, FDgrad, FDHess, kmax, n(i));
    Mfsol_ba(i,:) = fsol_ba;
    Mtsol_ba(i, :) =  tsol_ba;
    
    f_name = "fun_f_Banded";
    [fn_ba, ~, tn_ba] = test_NelderMead_random_points(f_name, n(i), kmax);
    Mfn_ba(i,:) = fn_ba;
    Mtn_ba(i, :) =  tn_ba;
end


fig2 = figure()
ymax = max(max(Mtn_ba))+ max(max(Mtn_ba))/100
ylim([0, ymax])

for j=1:3
   figure()
   hold on
   
   x = [0, n];
   vet = Mtsol_ba(:,j);
   y = [0, vet'];
   stairs(x, y)
   x = [0, n];
   vet = Mtn_ba(:,j);
   y = [0, vet'];
   stairs(x, y)

   legend("Exact Newton Method", "Nelder Mead Method")
   title('Banded trigonometric problem')
   xlabel('Dimension')
   ylabel('Computational time (s)')
   hold off
end
%}

%{
for i=1: length(n)
    f_Broyden = @fun_f_Broyden;
    f_gradient_Broyden = @grad_f_Broyden;
    f_Hessian_Broyden = @Hess_f_Broyden;

    [fsol_br, ~, ksol_br, tsol_br] = test_solver_random_points( ...
                            f_Broyden, f_gradient_Broyden, f_Hessian_Broyden, FDgrad, FDHess, kmax, n(i))
    Mfsol_br(i,:) = fsol_br;
    Mtsol_br(i, :) =  tsol_br;

    f_name = "fun_f_Broyden";
    [fn_br, kn_br, tn_br] = test_NelderMead_random_points(f_name, n(i), kmax)
    Mfn_br(i,:) = fn_br;
    Mtn_br(i, :) =  tn_br;

end


for j=1:2
   figure()
   hold on
   
   x = [0, n];
   vet = Mtsol_br(:,j);
   y = [0, vet'];
   stairs(x, y,'MarkerFaceColor')
   x = [0, n];
   vet = Mtn_br(:,j);
   y = [0, vet'];
   stairs(x, y, 'MarkerFaceColor')

   legend("Exact Newton Method", "Nelder Mead Method")
   title('Banded trigonometric problem')
   xlabel('Dimension')
   ylabel('Computational time (s)')
   hold off
end
%}
                           
