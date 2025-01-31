function [grad] = grad_f_Broyden(x)
    
    fk = @(x, k) (3-2*x(k))*x(k)-x(k-1)-2*x(k+1)+1;
    
    grad = zeros(length(x), 1);
    grad(1) = ((3-2*x(1))*x(1)-2*x(2)+1)*(3-4*x(1)) + fk(x, 2)*(-1);
    grad(2) = fk(x, 2)*(3-4*x(2)) + fk(x, 3)*(-1) + ((3-2*x(1))*x(1)-2*x(2)+1)*(-2);
    
    for i=3:length(x)-2
        grad(i) = fk(x, i)*(3-4*x(i)) + fk(x, i+1)*(-1) + fk(x, i-1)*(-2);
    end
    i = i+1;
    grad(i) = fk(x, i)*(3-4*x(i)) + ((3-2*x(i+1))*x(i+1)-x(i)+1)*(-1) + fk(x, i-1)*(-2);
    i = i+1;
    grad(i) = ((3-2*x(i))*x(i)-x(i-1)+1)*(3-4*x(i)) + fk(x, i-1)*(-2);
end