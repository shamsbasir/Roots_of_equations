function [x,numIts]=newton0(func,x,maxIts)
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

% SHAMSULHAQ BASIR 19.09.2018

if nargin < 3
  maxIts=100;      % default value if maxIts is omitted
end

EPSILON = 5.0e-5;
for numIts=1:maxIts
    
  [value,derivative]= func(x);
  increment = -(value/derivative);
  
  x = x + increment;
  
  errorEstimate = abs(increment);

  disp(strcat(num2str(numIts), ' x=', num2str(x), ' error estimate=', ...
    num2str(errorEstimate)));

  if errorEstimate<EPSILON
    return;
  end
end
% if get here, the Newton iteration has failed!
error('newton: maximum number of iterations exceeded.')
