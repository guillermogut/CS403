clear all
close all
clc
addpath('..\matlab_utils')

%================ q vectors for hw4 ==================
%q = [0,pi/2,0,pi/6,pi/2,0];%1A
%q = [0,2*pi/3,0,pi/3,pi/2,0];%1B
%q = [0,pi/2,pi/2,pi/6,pi/2,0]%2A
q = [0,pi/3,pi/4,pi/3,pi/2,0]%2B

%================ q vectors for hw4 ==================

th1 = q(1);
th2 = q(2);
th3 = q(3);
th4 = q(4);
th5 = q(5);
th6 = q(6);

R01 = e2O([th1,0,0]);
R12 = e2O([0,th2,0]);
R23 = e2O([0,th3,0]);
R34 = e2O([0,0,th4]);
R45 = e2O([0,th5,0]);
R56 = e2O([0,0,th6]);


T01 = SE3(R01,[0;0;0]);
T12 = SE3(R12,[0;0;.15]);
T23 = SE3(R23,[.3;0;0]);
T34 = SE3(R34,[.15;0;0]);
T45 = SE3(R45,[.1;0;0]);
T56 = SE3(R56,[.07;0;0]);
T6E = SE3(eye(3),[.05;0;0]);

T02 = T01*T12;
T03 = T02*T23;
T04 = T03*T34;
T05 = T04*T45;
T06 = T05*T56;
T0E = T06*T6E;

figure
clf;
grid on
axis equal

xlabel('$x$','interpreter','latex','fontsize',20)
ylabel('$y$','interpreter','latex','fontsize',20)
zlabel('$z$','interpreter','latex','fontsize',20)
xlim([-1, 1])
ylim([-1, 1])
zlim([-.5, .5])

view(60, 30);

drawCoordinate3DScale(eye(3),[0,0,0],.1);%frame0
drawCoordinate3DScale(T01(1:3,1:3),T01(1:3,4),.1);%frame1
drawCoordinate3DScale(T02(1:3,1:3),T02(1:3,4),.1);%frame2
drawCoordinate3DScale(T03(1:3,1:3),T03(1:3,4),.1);%frame3
drawCoordinate3DScale(T04(1:3,1:3),T04(1:3,4),.1);%frame4
drawCoordinate3DScale(T05(1:3,1:3),T05(1:3,4),.1);%frame5
drawCoordinate3DScale(T06(1:3,1:3),T06(1:3,4),.1);%frame6
hold on;
drawLine3D([0;0;0],T02(1:3,4));
drawLine3D(T02(1:3,4),T03(1:3,4));
drawLine3D(T03(1:3,4),T04(1:3,4));
drawLine3D(T04(1:3,4),T05(1:3,4));
drawLine3D(T05(1:3,4),T06(1:3,4));
drawLine3D(T06(1:3,4),T0E(1:3,4));

%drawCoordinate3D(T0E(1:3,1:3),T0E(1:3,4));%frameE

function x = e2O(angles)% convert ZYX Euler angle to an orientation matrix
thZ = angles(1);
thY = angles(2);
thX = angles(3);
Z = [cos(thZ)   -sin(thZ)   0;
    sin(thZ)     cos(thZ)   0;
    0           0           1];
Y = [cos(thY)  0  sin(thY);
       0       1      0;
    -sin(thY)  0  cos(thY)];
X = [1        0         0;
    0     cos(thX) -sin(thX);
    0     sin(thX) cos(thX)];    
x = Z*Y*X;
end    

function se = SE3(R,P)%get SE(3) group with rotation R and translation V
se = [ R    P ;
      0 0 0 1];
end

function product = mult2SE3(se1,se2)%not sure why we need this one
product = se1*se2;
end

function inverse = SE3Inverse(m)%inverse transformation matrix 
R = m(1:3,1:3);
P = m(1:3,4);
Rt = transpose(R);
negRt = -1*Rt;
inverse = [ Rt negRt*P;
            0 0 0 1];   
end
