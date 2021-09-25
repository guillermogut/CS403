clear all
close all
clc

addpath('..\matlab_utils')

%ZYX = (0.3, 0.2, 0.5),
%(0.7, π,π/2), 
%( π/3, 0, 0)

%%%------1-----%%%%

%(0.3, 0.2, 0.5)
orientation1 = e2O([0.3, 0.2, 0.5]);

%(0.7, π,π/2)
orientation2= e2O([0.7, pi,pi/2]);

%( π/3, 0, 0)
orientation3= e2O([ pi/3, 0, 0]);


T01 = SE3(orientation1,[0.4; 0.8; 1.2]);
T12 = SE3(orientation2,[-0.4; 0.5; 1.0]);
T23 = SE3(orientation3,[0.5; -0.8; 1.2]);


T02 = T01*T12;
T03 = T02*T23;

T30 = SE3Inverse(T03);

figure('position', [50, 0, 1800, 1000])
num_step = 50;
for i = 1:num_step
    
    clf;
    grid on
    axis equal

    xlim([-3, 4])
    ylim([-3, 4])
    zlim([0, 4])

    view(60, 30);


    drawCoordinate3D(eye(3),[0,0,0]);%frame 0
    drawCoordinate3D(T01(1:3,1:3),T01(1:3,4));%frame 1
    drawCoordinate3D(T02(1:3,1:3),T02(1:3,4));%frame 2
    drawCoordinate3D(T03(1:3,1:3),T03(1:3,4));%frame 3

    t = 0.1*i;
    T3E = SE3(orientation3,[0.1*sin(pi/3*t)+0.05;0.3*cos(pi/3*t)+0.08;sin(pi/3*t)+0.5]);

    T0E = T03*T3E;

    drawCoordinate3D(T0E(1:3,1:3),T0E(1:3,4));%origin frame 3
    
    pause(0.01);
    
    
end


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


