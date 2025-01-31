function [Hess] = Hess_f_Chained_Rosenbrock(x)

diagonal = zeros(1,length(x));
up_diagonal = zeros(1,length(x));
low_diagonal = zeros(1,length(x));

diagonal(1) = 100*12*x(1)^2 - 100*4*x(2) + 2;
up_diagonal(2) = -100*4*x(1);                       % the first term will we ignored


for i=2:length(x)-1
    diagonal(i) = 100*2 + 100*12*x(i)^2 - 100*4*x(i+1) + 2;
    up_diagonal(i+1) = -100*4*x(i);                  
    low_diagonal(i-1) = -100*4*x(i-1);
end

i = i+1;
diagonal(i) = 100*2;
low_diagonal(i-1) = -100*4*x(i-1);                  %the laste term will be ignore

Hess=spdiags([up_diagonal' diagonal' low_diagonal' ],[1 0 -1], length(x), length(x));

end