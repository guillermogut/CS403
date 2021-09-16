%%
addpath('..\matlab_utils')
syms th1 th2 th3 F1 F2 t1 t2 
% 

%th1,2 and 3 are angles
%length of arms


J = jacobian([l1*cos(th1) + l2* cos(th1 + th2) + l3* cos(th1+th2+th3); l1*sin(th1) + l2* sin(th1 + th2) + l3* sin(th1+th2+th3)],[th1;th2;th3]);


%relation of torque to forces using Jacobian

%torque * dq - F * dx = 0

%torque * dq = F * dx

%torque = F * (dx/dq)

%torque = J^T * F


