function [x,itCount] = secant(func, a, b)
% [x,itCount] = secant(func, a, b) uses secant method to find the
% root of func between a and b to tolerance of 1.0e-10
% a=left end point of interval
% b=right end point of interval
% x is the approximate root found
% itCount is the number of iterations required.
% ITMAX is the maximum number of iteration
% EPSILON is the tolerance for convergence

% SHAMSULHAQ BASIR
% 16.09.2018


EPSILON = 1.0e-10;
ITMAX = 1000;
for itCount = 1:ITMAX  
  x = b- (b-a)*func(b)/(func(b)-func(a));          
 if (abs(func(x)) <= EPSILON)
     return ; %converged;
 else
     a = b;
     b = x; 
 end
end
 error('Secant failed with too many iterations!')
 
 %% ANSWERS EXERCISE 8
 
 %% 1) 
 
 % code is created and modified accordingly
 
 %% 2)
 
 % EPSILON = 1.0e-10; ITMAX = 1000; are used
 
 %% 3)
 
 %disp is commented to turn off
 
 %% 4)
 
 % code is written
 
 %% 5) 
 
 % [x,n]=secant(@f0,0,3) is executed. root is 1, iteration no: 1
 
 %% 6)
 
 % excel spreadsheet is created and appended
 
 %% 7) 
 
 % as seen in the above table, secant method is faster for f1,f2,f4.
 % slower in f5 and f3. for f3 in [-1,2]does not converge.
 % In case of Bisection, f5 is close to the root which is 1. 
 