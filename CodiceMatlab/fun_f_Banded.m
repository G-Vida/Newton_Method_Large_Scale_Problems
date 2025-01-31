function [sum] = fun_f_Banded(x)

sum = 0;

%x_0(index 1) and x_(n-1)(index n) are 0
sum = 1*(1-cos(x(1))-sin(x(2))) + sum;
for i=2:length(x)-1
    sum = i*(1-cos(x(i))+sin(x(i-1))-sin(x(i+1))) + sum;
end
i = i+1;
sum = i*(1-cos(x(i))+sin(x(i-1))) + sum;

end