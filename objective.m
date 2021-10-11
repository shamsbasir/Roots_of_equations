function [f,J,F]=objective(x)
% [f,J,F]=objective(x) computes the objective function
% for a nonlinear least-square search for the exponent 
% and amplitudes for a fit to a time sequence.
% the time sequence is generated from 
%   v=exp(-ex1*t)*(ex3*sin(ex2*t)+ex4*cos(ex2*t));
% where ex1=0.15; ex2=2; ex3=1; ex4=3;
% f=grad F (column vector)
% J=Jacobian (4x4 matrix)
% F=objective function (scalar)

% $Id: objective.m,v 1.5 2012/09/10 01:14:36 mike Exp $
% M. M. Sussman

% these are the exact target values for x(1)-x(4)
ex1=0.15;
ex2=2;
ex3=1;
ex4=3;

t=linspace(0,40,800);
% v is the vector of time-series data
v=exp(-ex1*t).*(ex3*sin(ex2*t)+ex4*cos(ex2*t));

% initialize column vector f and matrix J
f=zeros(4,1);
J=zeros(4,4);

% objective function
F = sum((v - exp(-x(1)*t).*(x(3).*sin(x(2)*t) + x(4).*cos(x(2)*t))).^2);

% gradient
f(1) = sum( 2*(v - exp(-x(1)*t).*(x(3)*sin(x(2)*t) + x(4)*cos(x(2)*t))).*exp(-x(1)*t).*(x(3).*sin(x(2)*t) + x(4).*cos(x(2)*t)).*t);
f(2) = sum(-2*(v - exp(-x(1)*t).*(x(3)*sin(x(2)*t) + x(4)*cos(x(2)*t))).*exp(-x(1)*t).*(x(3).*cos(x(2)*t).*t - x(4).*sin(x(2)*t).*t));
f(3) = sum(-2*(v - exp(-x(1)*t).*(x(3)*sin(x(2)*t) + x(4)*cos(x(2)*t))).*exp(-x(1)*t).*sin(x(2)*t));
f(4) = sum(-2*(v - exp(-x(1)*t).*(x(3)*sin(x(2)*t) + x(4)*cos(x(2)*t))).*exp(-x(1)*t).*cos(x(2)*t));

% Jacobian matrix
J(1,1) = sum(2*t.^2 .* exp(-x(1)*t).^2 .* (x(3).*sin(x(2)*t) + x(4).*cos(x(2)*t)).^2  - 2 *(v - exp(-x(1)*t).*(x(3).*sin(x(2)*t) + ...
      x(4).*cos(x(2)*t))).*t.^2 .* exp(-x(1)*t).*(x(3).*sin(x(2)*t) + x(4).*cos(x(2)*t)));

J(1,2) = sum(-2*exp(-x(1)*t).^2 .*(x(3).*cos(x(2)*t).*t - x(4).*sin(x(2)*t).*t).*t.*(x(3).*sin(x(2)*t) + x(4).*cos(x(2)*t)) ...
     + 2*(v - exp(-x(1)*t).*(x(3).*sin(x(2)*t) + x(4).*cos(x(2)*t))).*t.*exp(-x(1)*t) ...
     .* (x(3).*cos(x(2)*t).*t - x(4).*sin(x(2)*t).*t));

J(1,3) = sum(-2*exp(-x(1)*t).^2 .*sin(x(2)*t).*t.*(x(3).*sin(x(2)*t) + x(4).*cos(x(2)*t)) ...
     + 2*(v - exp(-x(1)*t).*(x(3).*sin(x(2)*t) + x(4).*cos(x(2)*t))).*t.*exp(-x(1)*t).*sin(x(2)*t));

J(1,4) = sum(-2*exp(-x(1)*t).^2 .*cos(x(2)*t).*t.*(x(3).*sin(x(2)*t) + x(4).*cos(x(2)*t)) ...
     + 2*(v - exp(-x(1)*t).*(x(3).*sin(x(2)*t) + x(4).*cos(x(2)*t))).*t.*exp(-x(1)*t).*cos(x(2)*t));

J(2,1) = sum(-2*exp(-x(1)*t).^2 .*(x(3).*cos(x(2)*t).*t - x(4).*sin(x(2)*t).*t).*t.*(x(3).*sin(x(2)*t) + x(4).*cos(x(2)*t)) ...
     + 2*(v - exp(-x(1)*t).*(x(3).*sin(x(2)*t) + x(4).*cos(x(2)*t))).*t.*exp(-x(1)*t) ...
     .* (x(3).*cos(x(2)*t).*t - x(4).*sin(x(2)*t).*t));

J(2,2) = sum(2*exp(-x(1)*t).^2 .*(x(3).*cos(x(2)*t).*t - x(4).*sin(x(2)*t).*t).^2  - 2 ...
     .* (v - exp(-x(1)*t).*(x(3).*sin(x(2)*t) + x(4).*cos(x(2)*t))).*exp(-x(1)*t) ...
     .* (-x(3).*sin(x(2)*t).*t.^2  - x(4).*cos(x(2)*t).*t.^2 ));

J(2,3) = sum(2*exp(-x(1)*t).^2 .*sin(x(2)*t).*(x(3).*cos(x(2)*t).*t - x(4).*sin(x(2)*t).*t) ...
     - 2*(v - exp(-x(1)*t).*(x(3).*sin(x(2)*t) + x(4).*cos(x(2)*t))).*t.*exp(-x(1)*t).*cos(x(2)*t));

J(2,4) = sum(2*exp(-x(1)*t).^2 .*cos(x(2)*t).*(x(3).*cos(x(2)*t).*t - x(4).*sin(x(2)*t).*t) ...
     + 2*(v - exp(-x(1)*t).*(x(3).*sin(x(2)*t) + x(4).*cos(x(2)*t))).*t.*exp(-x(1)*t).*sin(x(2)*t));

J(3,1) = sum(-2*exp(-x(1)*t).^2 .*sin(x(2)*t).*t.*(x(3).*sin(x(2)*t) + x(4).*cos(x(2)*t)) ...
     + 2*(v - exp(-x(1)*t).*(x(3).*sin(x(2)*t) + x(4).*cos(x(2)*t))).*t.*exp(-x(1)*t).*sin(x(2)*t));

J(3,2) = sum(2*exp(-x(1)*t).^2 .*sin(x(2)*t).*(x(3).*cos(x(2)*t).*t - x(4).*sin(x(2)*t).*t) ...
     - 2*(v - exp(-x(1)*t).*(x(3).*sin(x(2)*t) + x(4).*cos(x(2)*t))).*t.*exp(-x(1)*t).*cos(x(2)*t));

J(3,3) = sum(2*exp(-x(1)*t).^2 .*sin(x(2)*t).^2);

J(3,4) = sum(2*exp(-x(1)*t).^2 .*cos(x(2)*t).*sin(x(2)*t));

J(4,1) = sum(-2*exp(-x(1)*t).^2 .*cos(x(2)*t).*t.*(x(3).*sin(x(2)*t) + x(4).*cos(x(2)*t)) ...
     + 2*(v - exp(-x(1)*t).*(x(3).*sin(x(2)*t) + x(4).*cos(x(2)*t))).*t.*exp(-x(1)*t).*cos(x(2)*t));

J(4,2) = sum(2*exp(-x(1)*t).^2 .*cos(x(2)*t).*(x(3).*cos(x(2)*t).*t - x(4).*sin(x(2)*t).*t) ...
     + 2*(v - exp(-x(1)*t).*(x(3).*sin(x(2)*t) + x(4).*cos(x(2)*t))).*t.*exp(-x(1)*t).*sin(x(2)*t));

J(4,3) = sum(2*exp(-x(1)*t).^2 .*cos(x(2)*t).*sin(x(2)*t));

J(4,4) = sum(2*exp(-x(1)*t).^2 .*cos(x(2)*t).^2);
