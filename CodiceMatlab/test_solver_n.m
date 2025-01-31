function [fsol, gradsol, ksol, tsol] = test_solver_n(f, grad, Hess, type, FDgrad, FDHess)


%Fixed parameters
rho = 0.5;                          
h = sqrt(eps);                     
c = 10^(-4);                     
tolgrad = 10^(-12);     
btmax = 100;                    
pcg_maxit = 100; 
kmax = 100000;


%Variable parameters
n = [100, 10, 10];               %i

 
fsol = [0, 0, 0];
gradsol = zeros(3, 1);
ksol = [0, 0, 0];
tsol = [0, 0, 0];
  

% Varying n
for i = 1:3
    
    dim = n(i) 
        
    %Starting point
    if type ==1
        x0 = zeros(dim,1);
        for k=1:dim
            if(mod(k,2)==1)
                x0(k) = -1.2;
            else
                x0(k) = 1;
            end
        end 
    end

    if type == 2
        x0 = ones(dim, 1);
    end
    
    if type == 3
        x0 = -ones(dim, 1);
    end
      
    
    tStart= cputime;
    [~, fk, gradfk_norm, k] = newton_backtrack(x0, f, grad, Hess, kmax,tolgrad, c, rho, btmax, FDgrad, FDHess, h, pcg_maxit);
    tEnd = cputime - tStart;

    
   
    fsol (i) = fk;
    gradsol(i) = gradfk_norm;
    ksol(i) = k;
    tsol(i)= tEnd;
    
    i = i+1;

end

end