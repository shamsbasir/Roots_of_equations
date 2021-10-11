function [x,itCount] = bisect(func, a, b)
% [x,itCount] = bisect=(func, a, b) uses bisection to find a
% root of func between a and b to tolerance of 1.0e-10
% a=left end point of interval
% b=right end point of interval
% func(a) and func(b) should be of opposite signs
% x is the approximate root found
% itCount is the number of iterations required.

% SHAMSULHAQ BASIR 
% 15.09.2018

%% 
EPSILON = 1.0e-10;
fa = func(a); 
fb = func(b);
n = ceil((log((b-a)/EPSILON)/log(2))-1)+1; % n from 3 exercise.

% check for not giving the wrong root
 if(sign(func(a))*sign(func(b))> 0)
     error('no root in this interval');
 end
 
for itCount = 1:n  % filled in using the formula from Exercise 3
    
  x = (b+a)/2;
  fx = func(x);

  if ( fx == 0 )
    return;  % found the solution exactly!
  elseif ( abs( b - x ) < EPSILON )
    return;  % satisfied the convergence criterion
  end

  if ( sign(fa) * sign(fx) <= 0 )
    b = x;
    fb = fx;
  else
    a = x; 
    fa = fx;
  end

end
 error('bisect failed with too many iterations!')
 
 
 