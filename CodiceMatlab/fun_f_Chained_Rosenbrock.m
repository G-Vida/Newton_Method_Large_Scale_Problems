function [sum] = fun_f_Chained_Rosenbrock(x)

sum = 0;
for i=2:length(x)
    sum = 100*(x(i-1)^2 - x(i))^2+(x(i-1)-1)^2 + sum;
end

end


