function [grad] = grad_f_Chained_Rosenbrock(x)

grad = zeros(length(x),1);

grad(1) =100*4*x(1)*(x(1)^2 - x(2)) + 2*(x(1)-1);

%forse 3
for i=2:length(x)-1
    grad(i) = -100*2*(x(i-1)^2 - x(i)) + 100*4*x(i)*(x(i)^2 - x(i+1)) + 2*(x(i)-1);
end

i = i+1;

grad(i) = -100*2*(x(i-1)^2 - x(i));

end
