function [fsol, ksol, tsol] = test_NelderMead_random_points(f_name, n, kmax)

rng(11);
xr0_1 = randi([-4, 5], 1, n)*0.75;
xr0_1 = xr0_1';

rng(15);
xr0_2 = randi([-4, 5], 1, n)*0.75;
xr0_2 = xr0_2';

rng(23);
xr0_3 = randi([-4, 5], 1, n)*0.75;
xr0_3 = xr0_3';

xr0 = [xr0_1, xr0_2, xr0_3];

fsol = [0, 0, 0];
ksol = [0, 0, 0];
tsol = [0, 0, 0];


for i = 1:3
        
    tStart= cputime;
    [fk,k] = NelderMead(f_name, xr0(:, i), n, kmax);
    tEnd = cputime - tStart;
    
    fsol(i) = fk;
    ksol(i) = k;
    tsol(i)= tEnd;
    
    i = i+1;

end

end