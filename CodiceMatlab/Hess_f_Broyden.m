function [Hess] = Hess_f_Broyden(x)

    diagonal = zeros(1,length(x));
    up_diagonal = zeros(1,length(x));
    low_diagonal = zeros(1,length(x));

    
    diagonal(1) = 6+24*x(1)^2-36*x(1)+8*x(1+1);
    up_diagonal(2) = -9+4*x(1+1)+8*x(1);                              

    %Saranno sempre 2:
    upup_diagonal = 2*ones(1, length(x));
    lowlow_diagonal = 2*ones(1, length(x));
    
for i=2:length(x)-2
    diagonal(i) = 10+24*x(i)^2-36*x(i)+4*x(i-1)+8*x(i+1);
    up_diagonal(i+1) = -9+4*x(i+1)+8*x(i);
    low_diagonal(i-1) = -9+4*x(i)+8*x(i-1); 
end

i = i+1;

diagonal(i) = 10+24*x(i)^2-36*x(i)+4*x(i-1)+8*x(i+1);
up_diagonal(i+1) = -9+4*x(i+1)+8*x(i);
low_diagonal(i-1) = -9+4*x(i)+8*x(i-1);

i = i+1;

diagonal(i) = 9+24*x(i)^2-36*x(i)+4*x(i-1);
low_diagonal(i-1) = -9+4*x(i)+8*x(i-1);

Hess=spdiags([upup_diagonal' up_diagonal' diagonal' low_diagonal' lowlow_diagonal'], [2 1 0 -1 -2], length(x), length(x));
end