
diagonal = zeros(1,10);
up_diagonal = zeros(1,10);
low_diagonal = zeros(1,10);

for i=1:10
    diagonal(i) = i;
    up_diagonal(i) = i;
    low_diagonal(i) = i;
end

Hess=spdiags([up_diagonal' diagonal' low_diagonal' ],[1 0 -1], 10, 10);
A = full(Hess)