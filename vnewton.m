function [x,numIts]=vnewton(func,x,maxIts)
% [x,numIts]=newton(func,x,maxIts)
% func is a function handle with signature [y,yprime]=func(x)
% on input, x is the initial guess
% on output, x is the final solution
% EPSILON is the convergence criteria.
% maxIts is the number of maximum iteration considered
% maxIts is an optional argument
% the default value of maxIts is 100
% numIts is the output of the iteration it took to find the root
% Newton's method is used to find x so that func(x)=0

% SHAMSULHAQ BASIR 30.09.2018

if nargin < 3
  maxIts=100;      % default value if maxIts is omitted
end

EPSILON = 5.0e-5;

increment=1;  % this is an arbitrary value
for numIts=1:maxIts
    
  [value,derivative]= func(x);
 
    oldIncrement=increment;
    
  increment = -(derivative\value);
  x = x + increment;
  r1 = norm(increment)/norm(oldIncrement) ;
  errorEstimate = norm(increment,2);

    if errorEstimate < EPSILON*(1-r1)
    return;
    end
    
end
% if get here, the Newton iteration has failed!
error('newton: maximum number of iterations exceeded.')
