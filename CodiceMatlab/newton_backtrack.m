function [xk, fk, gradfk_norm, k, non_positive, xseq, rel_diff, tot_time_precond, tot_time_nonprecond, exact_time] = ...
    newton_backtrack(x0, f, gradf , Hessf, kmax, tolgrad, c, rho, btmax, FDgrad, FDHess, h, pcg_maxit)
%
% [xk, fk, gradfk_norm, k, non_positive, xseq, rel_diff, tot_time_precond, tot_time_nonprecond] = 
%    newton_backtrack(x0, fk, gradf , Hessf, kmax, tolgrad, c, rho, btmax, FDgrad, FDHess, h, pcg_maxit)
%
% Function that performs the Newton optimization method with backtracking strategy.
%
% INPUTS:
% x0 = n-dimensional column vector;
% f = function handle that describes a function R^n->R;
% gradf = function handle that describes the gradient of f;
% Hessf = function handle that describes the Hessian of f;
% kmax = maximum number of iterations permitted;
% tolgrad = value used as stopping criterion w.r.t norm of the gradient;
% c = the factor of the Armijo condition that must be a scalar in (0,1);
% rho = fixed factor, lesser than 1, used for reducing alpha0;
% btmax = maximum number of steps for updating alpha during the backtracking strategy;
% FDgrad = 'fw' (FD Forward approx. for gradf),'c' (FD Centered approx. for gradf), any other string (usage of gradf)
% FDHess = 'fw' (FD Forward approx. for Hessf), 'Jfw' (Jacobian FD Forward approx. of Hessf), 'Jc' (Jacobian FD Centered approx. of Hessf), 
%    'MF' (Matrix Free implementation for solving Hessf(xk)pk=-gradf(xk)) any other string (usage of Hessf);
% h = approximation step for FD (if used);
% pcg_maxit = maximum number of iterations for the pcg solver.
%
% OUTPUTS:
% xk = the last x computed by the function;
% fk = the value f(xk);
% gradfk_norm = value of the norm of gradf(xk)
% k = index of the last iteration performed
% non_positive = a flag indicating non positivness of the Hessian of f (1 if non posive, 0 if positive)
% x_seq = sequence of x 
% rel_diff = relative difference (in infinity norm) of the fector founded with preconditioned and non preconditioned pcg
%tot_time_precond = computational time of finding the preconditioner + preconditioned pcg
%tot_time_nonprecond = computational time of non-preconditioned pcg
%exact_time = computational time of the direct exact method


switch FDgrad
    case 'fw'
        % OVERWRITE gradf WITH A F. HANDLE THAT USES findiff_grad
        % (with option 'fw')
        gradf = @(x) findiff_grad(f, x, norm(x)*h, 'fw');
        
    case 'c'
        % OVERWRITE gradf WITH A F. HANDLE THAT USES findiff_grad
        % (with option 'c')        
        gradf = @(x) findiff_grad(f, x, norm(x)*h, 'c');
        
    otherwise
        % WE USE THE INPUT FUNCTION HANDLE gradf  
end


%CASO DIMESIONE ECCESSIVA


%IN CASE OF APPROXIMATED GRADIENT, IT IS BETTER TO NOT APPROXIMATE Hessf WITH THE JACOBIAN!
if isequal(FDgrad, 'fw') || isequal(FDgrad, 'c')
    switch FDHess
        case 'fw'
            % WRITE/OVERWRITE Hessf WITH A F. HANDLE THAT USES findiff_Hess
            Hessf = @(x) findiff_Hess(f, x, norm(x)*h); 
        case 'MF'
            % DEFINE a f. handle for the product of Hessf * p USING THE GRADIENT
            Hessf_pk = @(x, p) (gradf(x + (norm(x)*h) * p) - gradf(x)) / (norm(x)*h);
        otherwise
            % WE USE THE INPUT FUNCTION HANDLE Hessf

    end
else
    switch FDHess
        case 'fw'
            % OVERWRITE Hessf WITH A F. HANDLE THAT USES findiff_Hess
            Hessf = @(x) findiff_Hess(f, x, norm(x)*h);                    %h variabile
        
        case 'Jfw'
            % OVERWRITE Hessf WITH A F. HANDLE THAT USES findiff_J
            % (with option 'fw')
            Hessf = @(x) findiff_J(gradf, x, norm(x)*h, 'fw');

        case 'Jc'
            % OVERWRITE Hessf WITH A F. HANDLE THAT USES findiff_J
            % (with option 'c')
            Hessf = @(x) findiff_J(gradf, x, norm(x)*h, 'c');

        case 'MF'
            % DEFINE a f. handle for the product of Hessf * p USING THE
            % GRADIENT
            Hessf_pk = @(x, p) (gradf(x + (norm(x)*h) * p) - gradf(x)) / (norm(x)*h);
        otherwise
            % WE USE THE INPUT FUNCTION HANDLE Hessf
    end
end



% Function handle for the armijo condition
farmijo = @(fk, alpha, gradfk, pk) fk + c * alpha * gradfk' * pk;

% Initializations
xseq = zeros(length(x0), kmax);
rel_diff = zeros(1,kmax);
ichol_time = zeros(1, kmax);
pcg_time = zeros(1, kmax);
tot_time_precond = zeros(1, kmax);
tot_time_nonprecond = zeros(1, kmax);
exact_time = zeros(1, kmax);

xk = x0;
fk = f(xk);
k = 0;
gradfk = gradf(xk);
gradfk_norm = norm(gradfk);
non_positive = 0;

tStart = tic;

while k < kmax && gradfk_norm >= tolgrad
    %Compute the descent direction as solution of Hessf(xk) p = - gradf(xk)
    switch FDHess
        case 'MF'
            % ITERATIVE METHOD: pcg with a function handle
            Hessfk_pk = @(p) Hessf_pk(xk, p);
            
            pk = pcg(Hessfk_pk, -gradfk, tolgrad, pcg_maxit);
            
        otherwise       
            % Check if the Hessian matrix is symmetric positive definitive:
            A = Hessf(xk);
            tic
            [~, FLAG] = chol(A);
            toc
            
            %If the Hesssian is not pos.def., we preprocess it adding a correction
            if FLAG ~= 0   
                non_positive = 1;
                Bk = A;
                I = spdiags([ones(length(xk), 1)],  [ 0 ], length(xk), length(xk));
                beta = norm(A,"fro");
                tauk = 0;
                
                while FLAG ~= 0    
                    tauk = max(2*tauk, beta/2);
                    Bk = A + tauk*I;
                    [~, FLAG] = chol(Bk);
                end
                
                Bk = sparse(Bk);
                
                if length(x0) > 10^1
                    
                    disp("BGULOOOO")
                    try
                        tic
                        R = ichol(Bk);
                        %toc
                        ichol_time(k+1) = toc;

                        tic
                        %We use R and R' as preconditioners for pcg
                        [pk, ~, ~, m] = pcg(Bk, -gradfk, tolgrad, pcg_maxit, R, R');
                        %toc
                        pcg_time(k+1) = toc;

                    catch
                        q = max(sum(abs(Bk),2)./diag(Bk))-2;
                        tic
                        R = ichol(Bk, struct('type','ict','droptol',1e-3,'diagcomp',q));
                        %toc
                        ichol_time(k+1) = toc;

                        tic
                        [pk, ~, ~, m] = pcg(Bk, -gradfk, tolgrad, pcg_maxit, R, R');    
                        %toc
                        pcg_time(k+1) = toc;
                    end
                    tot_time_precond(k+1) = ichol_time(k+1) + pcg_time(k+1);

                    % COMPARISON--------------------------------
                    tic
                    [tk, ~, ~, n] = pcg(Bk, -gradfk, tolgrad, pcg_maxit);
                    %toc
                    tot_time_nonprecond(k+1) = toc;

                    rel_diff(k+1) =  norm(pk-tk, Inf)/norm(pk, Inf);

                else
                    disp("EXACT")
                    tic
                    pk = Bk\gradfk;
                    exact_time(k+1) = toc;
                end
            else
                % We use R and R' as preconditioners for pcg if ichol a negative pivot will 
                % be found the function brake down and we use shifted incomplete cholesky
                % that guarantees R*R' is a diagonal dominant matric(positive definite)

                if length(x0) > 10^1
                    A = sparse(A);
                    try
                        tic
                        R = ichol(A);
                        %toc
                        ichol_time(k+1) = toc;

                        tic
                        [pk, ~, ~, m] = pcg(A, -gradfk, tolgrad, pcg_maxit, R, R');
                        %toc
                        pcg_time(k+1) = toc;

                    catch
                        q = max(sum(abs(A),2)./diag(A));
                        tic
                        R = ichol(A, struct('type','ict','droptol',1e-3,'diagcomp',q));
                        %toc
                        ichol_time(k+1) = toc;
                        
                        tic
                        [pk, ~, ~, m] = pcg(A, -gradfk, tolgrad, pcg_maxit, R, R');    %Solving with A (original Hessf)
                        %toc
                        pcg_time(k+1) = toc;
                    end
                    tot_time_precond(k+1) = ichol_time(k+1) + pcg_time(k+1);
                    

                    % COMPARISON--------------------------------
                    tic
                    [tk, ~, ~, n] = pcg(A, -gradfk, tolgrad, pcg_maxit);
                    toc
                    tot_time_nonprecond(k+1) = toc;
                    
                    pcg_iterations(k+1) = n;
                    pcg_precond_iterations(k+1) = m;                
                    rel_diff(k+1) =  norm(pk-tk, Inf)/norm(pk, Inf);
                  
                else
                    disp("EXACT")
                    tic
                    pk = A\-gradfk;  %Solving with A (original Hessf)
                    %toc
                    pcg_time(k+1) = toc;
                end
            end
    tEnd = toc(tStart)
    end


    
    
    % Reset the value of alpha
    alpha = 1;

    % Compute the candidate new xk
    xnew = xk + alpha * pk;
    % Compute the value of f in the candidate new xk
    fnew = f(xnew);
    bt = 0;
    % Backtracking strategy: 
    % 2nd condition is the Armijo condition not satisfied

    while bt < btmax && fnew > farmijo(fk, alpha, gradfk, pk)
        %Reduce the value of alpha
        alpha = rho * alpha;
        %Update xnew and fnew w.r.t. the reduced alpha
        xnew = xk + alpha * pk;
        fnew = f(xnew);
        %Increase the counter by one
        bt = bt + 1;
    end
    
    % Update xk, fk, gradfk_norm
    xk = xnew;
    fk = fnew;
    gradfk = gradf(xk);
    gradfk_norm = norm(gradfk);

    
    % Increase the step by one
    k = k + 1;

    % Store current xk in xseq
    xseq(:, k) = xk;
end
% "Cut" xseq
xseq = xseq(:, 1:k);
end