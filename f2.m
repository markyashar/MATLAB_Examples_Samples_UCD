function y = f2(x,par)
% this function evalueas a quadratic function
% This function is used in homework problem #1 in 
% homework assignment #1.
y = par.a*x.^2 + par.b*x + par.c;
% the "." in "x.^2" is crucial.  it tells mat lab to produce a new vector
% the same length as x with the ith entry equal to x(i)^2. Without the "."
% matalb would produce an error.  We could also write x.*x to get the same
% result. x*x would genereate errors but x*x' would give the dot product
% between x and iteslf. 