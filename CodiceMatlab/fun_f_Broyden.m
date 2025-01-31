function [sum] = fun_f_Broyden(x)
    sum = 1/2*(((3-2*x(1))*x(1)-2*x(2)+1)^2);
    
    for i=2:length(x)-1
        sum = 1/2*(((3-2*x(i))*x(i)-x(i-1)-2*x(i+1)+1)^2)+sum;
    end
    i = i+1;
    sum = 1/2*(((3-2*x(i))*x(i)-x(i-1)+1)^2)+sum;
end
