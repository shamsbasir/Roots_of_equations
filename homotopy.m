function [f,J,F]=homotopy(x,p,x0)
% [f,J,F]=homotopy(x,p,x0)
% computes the homotopy or Davidenko objective function
% for 0<=p<=1

[f1,J1,F1]=objective(x);
[f2,J2,F2]=easy_objective(x,x0);
f=p*f1+(1-p)*f2;
J=p*J1+(1-p)*J2;
F=p*F1+(1-p)*F2;