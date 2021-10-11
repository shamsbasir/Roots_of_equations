function [f,J,F]=easy_objective(x,x0)
% [f,J,F]=easy_objective(x,x0)
% f=grad F (column vector)
% J=Jacobian (4x4 matrix)
% F=objective function (scalar)


% shamsulhaq basir     '07-Oct-2018'

if norm(size(x)-size(x0)) ~= 0
  error('easy_objective: x and x0 must be compatible.')
end
F=sum((x-x0).^2);


% f(k)=derivative of F with respect to x(k)

f=zeros(4,1);
f(1,1)=sum(2*(x(1)-x0(1)));

f(2,1)=sum(2*(x(2)-x0(2)));

f(3,1)=sum(2*(x(3)-x0(3)));

f(4,1)=sum(2*(x(4)-x0(4)));


%J(k,ell)=derivative of f(k) with respect to x(ell)
J=diag([2,2,2,2]);