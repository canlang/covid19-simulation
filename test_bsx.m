fun = @(a,b) linspace(a,b,10);
% A = bsxfun(fun, 1:3,4:6);
A = arrayfun(fun, 1:3,4:6,'UniformOutput',false);

