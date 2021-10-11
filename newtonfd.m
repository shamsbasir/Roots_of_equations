function [x,numIts]=newtonfd(f,x,dx,maxIts)
% [x,numIts]=newtonfd(f,x,dx,maxIts)
% f is a function handle with signature [y,yprime]=f(x)
% on input, x is the initial guess
% on output, x is the final solution
% dx is step size. 
% EPSILON is the convergence criteria.
% maxIts is the number of maximum iteration considered
% maxIts is an optional argument
% the default value of maxIts is 100
% numIts is the output of the iteration it took to find the root
% Newton's method is used to find x so that f(x)=0

% SHAMSULHAQ BASIR 19.09.2018

if nargin < 4
  maxIts=100;                       % default value if maxIts is omitted
end

EPSILON = 5.0e-5;

increment=1;                        % this is an arbitrary value
for numIts=1:maxIts
    
  [value]= f(x);
 derivative = (f(x+dx)-f(x))/dx ;   % evaluation of derivative
 oldIncrement=increment;
    
  increment = -(value/derivative);
  x = x + increment;
  r1 = abs(increment)/abs(oldIncrement) ;
  r2 = abs(increment)/abs(oldIncrement)^2 ;
  errorEstimate = abs(increment);

%     disp(strcat(num2str(numIts), ' x=', num2str(x), ' error estimate=', ...
%       num2str(errorEstimate)));
% 
%    disp(strcat(num2str(numIts), ' r1=', num2str(r1), ' r2=', ...
%     num2str(r2)));

    if errorEstimate < EPSILON*(1-r1)
    return;
    end
    
end
% if get here, the Newton iteration has failed!
error('newton: maximum number of iterations exceeded.')

