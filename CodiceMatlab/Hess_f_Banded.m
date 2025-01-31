function [Hess] = Hess_f_Banded(x)

diagonal = zeros(1,length(x));

%x_0(index 1) and x_(n-1)(index n) are 0
diagonal(1) = 1*(cos(x(1)) - 2*sin(x(1)));
for i=2:length(x)-1
    diagonal(i) = i*cos(x(i)) + (i-1)*sin(x(i)) - (i+1)*sin(x(i));
end
i = i+1;
diagonal(i) = i*cos(x(i)) + (i-1)*sin(x(i));

Hess = spdiags([ diagonal' ],  [ 0 ], length(x), length(x));
end