function [x,numIts]=newton_sqrt(a,maxIts)
% calculates square root of the number a;
% maxIts is the maximum iteration number;
% EPSILON  is the convergence criteria;

% SHAMSULHAQ BASIR
% 25.09.2018

if nargin < 2
  maxIts=100;      % default value if maxIts is omitted
end

EPSILON = 10^-15;
x = 0.5*(a+1);
for numIts=1:maxIts
    old_x = x;
    x = 0.5*(x+a/x);
    trueError = abs(x-old_x);
    errorEstimate= abs(trueError/x);
    if errorEstimate < EPSILON
    return;
    end
end
% if get here, the Newton_sqrt iteration has failed!
error('newton_sqrt: maximum number of iterations exceeded.')
