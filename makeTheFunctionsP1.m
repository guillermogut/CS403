close all
clc

%code for questions 1 and 2
syms th1 th2 dth1 dth2 ddth1 ddth2 real
syms m1 m2 c1 c2 l1 l2 I1 I2 tau1 tau2 g real


q = [th1; th2];
dq = [dth1; dth2];
ddq = [ddth1; ddth2];
u = [tau1; tau2];


% params
p = [m1; m2; I1; I2; l1; l2; c1; c2; g];


% Function to take derivatives
ddt = @(r) jacobian(r, [q; dq])*[dq;ddq];

rCM1 = [c1*sin(th1);-c1*cos(th1)];
rB = [l1*sin(th1);-l1*cos(th1)];
rCM2 = [l1*sin(th1) + c2*sin(th1 +th2);
        -l1*cos(th1) - c2*cos(th1 +th2)];
rC = [l1*sin(th1) + l2*sin(th1 +th2);
        -l1*cos(th1) - l2*cos(th1 +th2)];

%derivatives
drCM1 = ddt(rCM1);
drB = ddt(rB);
drCM2 = ddt(rCM2);
drC = ddt(rC);

% kinetic Energy
T1 = (1/2)*m1*dot(drCM1, drCM1) + (1/2) * I1 * (dth1)^2;
T2 = (1/2)*m2*dot(drCM2, drCM2) + (1/2) * I2 * (dth1 + dth2)^2;

% gravity potential energy
V1 = m1*g*dot(rCM1, -([0;-1]));
V2 = m2*g*dot(rCM2, -([0;-1]));

% lagrangian
L = T1+T2-V1-V2;

M2Q = @(M,w) simplify(jacobian(w,dq)'*(M));

% tau1 and tau2 to gen. forces
Qt1 = M2Q(tau1*[0;0;1], dth1*[0;0;1]);
Qt2 = M2Q(tau2*[0;0;1], (dth1+dth2)*[0;0;1]);
Qt3 = M2Q(tau2*[0;0;1], (-dth1)*[0;0;1]);
Q = Qt1 + Qt2 + Qt3;

% implicit 
g = ddt(jacobian(L,dq)') - jacobian(L,q)' - Q;

% get ddq
A = jacobian(g,ddq);
b = A*ddq - g;
ddq2 = A\b;

% z and dz
z = [q; dq];
dz = [dq; ddq2];

% Total Energy
E = T1 + T2 + V1 + V2;

keypoints = [rB; rC];

matlabFunction(A, 'file', ['calc_A'], 'vars', {z p});
matlabFunction(b, 'file', ['calc_b'], 'vars', {z u p});
matlabFunction(keypoints, 'file', ['keypoints'], 'vars', {z p});
matlabFunction(E, 'file', ['energy'], 'vars', {z u p});




      