function [x,numIts]=snewton2(func,x,maxIts)
% function [x,numIts]=snewton1(func,x,maxIts)
% func is a function handle with signature [y,yprime]=func(x)
% on input, x is the initial guess
% on output, x is the final solution
% EPSILON is the convergence criteria.
% maxIts is the number of maximum iteration considered
% maxIts is an optional argument
% the default value of maxIts is 100
% numIts is the output of the iteration it took to find the root
% Newton's method is used to find x so that func(x)=0
% alpha  softening factor = 1/(1+10*abs(increment))
% SHAMSULHAQ BASIR 30.09.2018

if nargin < 3
  maxIts=100;      % default value if maxIts is omitted
end

EPSILON = 5.0e-5;
increment=1;  % this is an arbitrary value
for numIts=1:maxIts
    
  [value,derivative]= func(x);
 
    oldIncrement=increment;
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % calculation of softening factor 
  alpha = 1;
  increment = (derivative\value);
    for (j=1:10)
      
        xl = x + -alpha*increment;
  
        if(norm(func(xl))<= norm(func(x)) || alpha < 1/1024)
      
                x = xl;
             else
                     alpha = alpha/2;
            end
  
     end
  %%%%%%%%%%%%%%%%%%%%%%%
  
  r1 = norm(increment,2)/norm(oldIncrement,2) ;
  errorEstimate = norm(increment,2);

    if errorEstimate < EPSILON*(1-r1)
    return;
    end
    
end
% if get here, the Newton iteration has failed!
error('newton: maximum number of iterations exceeded.')
