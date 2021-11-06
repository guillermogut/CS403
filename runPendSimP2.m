close all
clear all
clc
%code for questions 3 and 4
tau1 = 0;
tau2 = 0;
th1 = 3;
th2 = 0;
dth1 = 0;
dth2 = 0;
m1 = 1;
m2 = 1;
I1 = 0.05;
I2 = 0.05;
l1 = 1;
l2 = 0.5;
c1 = 0.5;
c2 = 0.25;
k1 = 10;
k2 = 10;
k3 = 50;
th1_0 = 0;
th2_0 = 0;
l_0 = 0;

g = 9.8;
q = [th1; th2];
dq = [dth1; dth2];
u = [tau1; tau2];

r0 = [0; 0.5];

z = [q; dq];

p = [m1; m2; I1; I2; l1; l2; c1; c2; k1; k2; k3; th1_0; th2_0; l_0; r0; g];

th1list = [];
th2list = [];
totEnergy = [];

for i = 0:700
   
    dt = 0.01;

    A = calc_AS(z, p);
    b = calc_bS(z, u, p);
    ddq2 = A\b;

    total_energy = totEnergyS(z, [-z(3); -z(4)], p);
    totEnergy = [totEnergy;  total_energy];
    
    z = z + dt*[z(3:4) + dt*ddq2; ddq2];
   
    pos = keypoints(z, p);
    rB = pos(1:2);
    rC = pos(3:4);

    th1list = [th1list;z(1)];
    th2list = [th2list; z(2)];
    

    clf

    hold on
    drawLine2D([0; 0], rB);
    drawLine2D(rB, rC);
    drawLine2D(r0, rC);
    
    xlabel('$i\  (m)$','interpreter','latex','fontsize',15)
    ylabel('$j\  (m)$','interpreter','latex','fontsize',15)
    axis equal

    axis([-2 2 -2 2])
    grid on
    hold off

    pause(0.0001);
end

figure
plot(th1list);
grid on
xlim([0,700]);
xlabel('$t\  (s)$','interpreter','latex','fontsize',15)
ylabel('$\th1\  (rad)$','interpreter','latex','fontsize',15)

figure
plot(th2list);
grid on
xlim([0,700]);
xlabel('$t\  (s)$','interpreter','latex','fontsize',15)
ylabel('$\th2\  (rad)$','interpreter','latex','fontsize',15)

figure
plot(totEnergy);
grid on
xlim([0,700]);
xlabel('$t\  (s)$','interpreter','latex','fontsize',15)
ylabel('$Total\ Energy\  (J)$','interpreter','latex','fontsize',15)