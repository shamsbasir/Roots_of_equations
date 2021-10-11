function [y,yprime] = f3 ( x )
% Exercise 1: 
% y = x*exp(-x)computes the function value at x;
% yprime = exp(-x) - x*exp(-x) is the derivative of y;

% SHAMSULHAQ BASIR
% 19.09.2018
y = x*exp(-x);
yprime = exp(-x) - x*exp(-x);
end



