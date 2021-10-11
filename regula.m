function [x,itCount] = regula(func, a, b)
% [x,itCount] = regula(func, a, b) uses secant method to find the
% root of func between a and b to tolerance of 1.0e-10
% a=left end point of interval
% b=right end point of interval
% x is the approximate root found
% itCount is the number of iterations required.
% ITMAX is the maximum number of iteration
% EPSILON is the tolerance for convergence

% SHAMSULHAQ BASIR
% 16.09.2018

EPSILON = 1.0e-30;
ITMAX = 1000;
for itCount = 1:ITMAX  
  x = b- (b-a)*func(b)/(func(b)-func(a));          
 if (abs(func(x))<= EPSILON)
     return ; %converged;
 elseif(sign(func(a))*sign(func(b)) >=0 )
     a = b; 
 end
  b = x;
end
 error('regula failed with too many iterations!')
 