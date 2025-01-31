function [grad] = grad_f_Banded(x)

grad = zeros(1,length(x));

%x_0(index 1) and x_(n-1)(index n) are 0
grad(1) = 1*(sin(x(1)) + 2*cos(x(1)));
for i=2:length(x)-1
    grad(i) = i*sin(x(i)) - (i-1)*cos(x(i)) + (i+1)*cos(x(i));
end
i = i+1;
grad(i) = i*sin(x(i)) - (i-1)*cos(x(i));
grad = grad';
end