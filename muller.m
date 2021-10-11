function [result, itCount] = muller( func, a, b)
%[x,itCount] = muller(func, a, b) uses secant method to find the
% root of func between a and b to tolerance of 1.0e-10
% a=left end point of interval
% b=right end point of interval
% x is the approximate root found
% itCount is the number of iterations required.
% ITMAX is the maximum number of iteration
% EPSILON is the tolerance for convergence

%  Muller's method goes one better and approximates the function 
%  by the quadratic function passing through the two endpoints and 
%  the current iterate. The next iterate is then taken to be the root of 
%  this quadratic that is closest to the current iterate. 

% SHAMSULHAQ BASIR
% 16.09.2018
  EPSILON = 1.0e-10;
  ITMAX = 100;
  % start of the iteration
  x0 = a;
  x2 = b;
  x1 = .51*x0 + .49*x2;
  for itCount = 1:ITMAX 
% calculating ordinates of the points
  y0 = func(x0);
  y1 = func(x1);
  y2 = func(x2);
% y(x) =A(x-x2)^2 +B(x-x2)^2 +C approximation of the method
  A = ( (y0 - y2) * (x1 - x2) - (y1 - y2) * (x0 - x2) ) / ...
        ( (x0 - x2) * ( x1- x2) * (x0 - x1) );
  B = ( (y1 - y2) * (x0 - x2)^2 - (y0 - y2) * (x1 - x2)^2 ) / ...
        ( (x0 - x2) * (x1 - x2) * (x0 - x1) );
  C = y2;
  
%  If the polynomial has real roots, find them and choose the one 
%  closer to x2, otherwise choose something reasonable:
  
   if A ~= 0

    disc = B*B - 4.0*A*C;
    disc = max( disc, 0.0 );
         
    q1 = (B + sqrt(disc) );
    q2 = (B - sqrt(disc) );

    if abs(q1) < abs(q2) 
      dx = -2.0*C/q2;
    else
      dx = -2.0*C/q1;
    end

      elseif B ~= 0
        dx = -C/B;
            else
        error(['muller: algorithm broke down at itCount=', num2str(itCount)])
        
   end
  % Discarding the point (x0, y0) and add the new point   
  
  x0 = x1;
  y0 = y1;

  x1 = x2;
  y1 = y2;

  x2 = x1 + dx;
  y2 = func(x2);
  
  % checking if convergence criteria is reached
  
  if (abs(y2) <= EPSILON)
      result = x2;
      return;
  end
  end
error('muller failed with too many iterations!')  
  
  
  %%  ANSWERS TO EXERCISE 10 EXTRA-CREDIT
  
  %% 1) 
  
  % code is developed 
  
  %% 2) 
  
  % [x,n]=muller(@f0,0,2) is executed. root is reached in 1 iteration
  
  %% 3) 
  
  % [x,n]=muller(@f1,0,5) is executed. root is 3, reached in 1 iteration
  
  %% 4)
  
  % table is appended and it can be seen the method is very powerful and
  % fast for f3, muller diverges because of the fact interval at each
  % iteration need to be a change of sign to have a root. 
  
  
end

