clear all
close all
clc
addpath('..\matlab_utils')



q = [0,-pi/6,0,0,0,0];%starting position

th1 = q(1);
th2 = q(2);
th3 = q(3);
th4 = q(4);
th5 = q(5);
th6 = q(6);
th =[th1; th2; th3;th4;th5;th6];

S = {};

S{1} = zeros(6,1); S{1}(3) = 1;
S{2} = zeros(6,1); S{2}(2) = 1;
S{3} = zeros(6,1); S{3}(2) = 1;
S{4} = zeros(6,1); S{4}(1) = 1;
S{5} = zeros(6,1); S{5}(2) = 1;
S{6} = zeros(6,1); S{6}(1) = 1;

P = {};

P{1} = [0;0;.15];
P{2} = [.3;0;0];
P{3} = [.15;0;0];
P{4} = [.1;0;0];
P{5} = [.07;0;0];
P{6} = [.05;0;0];


tdRot = [0 -1  0;
         1  0  0;
         0  0  1];

%destination for robot arm EE
T_d = SE3(tdRot,[.2;.31;.2]);

R01 = SO3([th1,0,0]);
R12 = SO3([0,th2,0]);
R23 = SO3([0,th3,0]);
R34 = SO3([0,0,th4]);
R45 = SO3([0,th5,0]);
R56 = SO3([0,0,th6]);


T01 = SE3(R01,[0;0;0]);
T12 = SE3(R12,[0;0;.15]);
T23 = SE3(R23,[.3;0;0]);
T34 = SE3(R34,[.15;0;0]);
T45 = SE3(R45,[.1;0;0]);
T56 = SE3(R56,[.07;0;0]);
T6E = SE3(eye(3),[.05;0;0]);

T0 = {};
T02 = T01*T12;
T03 = T02*T23;
T04 = T03*T34;
T05 = T04*T45;
T06 = T05*T56;
T0E = T06*T6E;

T0{1} = T01;
T0{2} = T02;
T0{3} = T03;
T0{4} = T04;
T0{5} = T05;
T0{6} = T06;
T0{7} = T0E;
J = zeros(6,6);

figure('position', [850, 300, 450, 350]);
%===============================================================================================================
steps = 100;
for time = 1:steps
    clf;
    view(60, 30);  
    T0 = makeDaTransformations(th,P);
    
    J = jacobian(T0);
    %calculate the error
    
    error = zeros(6,1);
    error(1:3) = SO3_to_so3_Vector(T_d(1:3,1:3)*T0{7}(1:3,1:3)');
    error(4:6) = T_d(1:3,4)- T0{7}(1:3,4);
    
    %get theta for the next step
    th = th + 0.1*pinv(J)*error;
    
    grid on
    hold on
    axis equal
    drawCoordinate3DScale(eye(3),[0,0,0],.1);%frame0
    drawCoordinate3DScale(T0{1}(1:3,1:3),T0{1}(1:3,4),.05);%frame1
    drawCoordinate3DScale(T0{2}(1:3,1:3),T0{2}(1:3,4),.05);%frame2
    drawCoordinate3DScale(T0{3}(1:3,1:3),T0{3}(1:3,4),.05);%frame3
    drawCoordinate3DScale(T0{4}(1:3,1:3),T0{4}(1:3,4),.05);%frame4
    drawCoordinate3DScale(T0{5}(1:3,1:3),T0{5}(1:3,4),.05);%frame5
    drawCoordinate3DScale(T0{6}(1:3,1:3),T0{6}(1:3,4),.05);%frame6
    drawCoordinate3DScale(T0{7}(1:3,1:3),T0{7}(1:3,4),.1)
    hold on;
    drawLine3D([0;0;0],T0{2}(1:3,4));
    drawLine3D(T0{2}(1:3,4),T0{3}(1:3,4));
    drawLine3D(T0{3}(1:3,4),T0{4}(1:3,4));
    drawLine3D(T0{4}(1:3,4),T0{5}(1:3,4));
    drawLine3D(T0{5}(1:3,4),T0{6}(1:3,4));
    drawLine3D(T0{6}(1:3,4),T0{7}(1:3,4));

    xlabel('$x$','interpreter','latex','fontsize',20)
    ylabel('$y$','interpreter','latex','fontsize',20)
    zlabel('$z$','interpreter','latex','fontsize',20)
    xlim([-.5,.5])
    ylim([-.5, 1])
    zlim([-.2, .5])

    pause(0.00001);
    nor = norm(error);
    if nor < .1
        hold on;
        break
    end
    hold off;
end
hold on;

function J = jacobian(T0)%make the good ol' jacobian
    %axese of rotation
    S = {};
    S{1} = zeros(6,1); S{1}(3) = 1;
    S{2} = zeros(6,1); S{2}(2) = 1;
    S{3} = zeros(6,1); S{3}(2) = 1;
    S{4} = zeros(6,1); S{4}(1) = 1;
    S{5} = zeros(6,1); S{5}(2) = 1;
    S{6} = zeros(6,1); S{6}(1) = 1;
    J = zeros(6,6);
    %build the jacobian matrix
    for x = 1:6
        
        J(:,x) = adjoint(SE3Inverse(SE3Inverse(T0{x})*T0{7})) * S{x};

    end

    %change it to global frame 
    J = [T0{7}(1:3,1:3) zeros(3,3);zeros(3,3) T0{7}(1:3,1:3)]*J;
end


function  T0 =  makeDaTransformations(thetas,ps)% makes all transformation matricese

th1 = thetas(1);
th2 = thetas(2);
th3 = thetas(3);
th4 = thetas(4);
th5 = thetas(5);
th6 = thetas(6);

R01 = SO3([th1,0,0]);
R12 = SO3([0,th2,0]);
R23 = SO3([0,th3,0]);
R34 = SO3([0,0,th4]);
R45 = SO3([0,th5,0]);
R56 = SO3([0,0,th6]);


T01 = SE3(R01,[0;0;0]);
T12 = SE3(R12,ps{1});
T23 = SE3(R23,ps{2});
T34 = SE3(R34,ps{3});
T45 = SE3(R45,ps{4});
T56 = SE3(R56,ps{5});
T6E = SE3(eye(3),ps{6});

T0 = {};
T02 = T01*T12;
T03 = T02*T23;
T04 = T03*T34;
T05 = T04*T45;
T06 = T05*T56;
T0E = T06*T6E;

T0{1} = T01;
T0{2} = T02;
T0{3} = T03;
T0{4} = T04;
T0{5} = T05;
T0{6} = T06;
T0{7} = T0E;


end
function x = SO3_to_so3_Vector(R)

    so3 = logm(R);
    vectSMOLs03 = [so3(3,2);so3(1,3);so3(2,1)];
    x = vectSMOLs03;
end

function x = adjoint(T)%make a transformation into an adjoint matrix, T_EE to i

x = zeros(6,6);
%px py pz
P = [T(1,4),T(2,4),T(3,4)];

PSkew = [0  -P(3)  P(2);
        P(3)  0   -P(1);
        P(2) P(1)   0 ];

%build the return matrix
x(1:3,1:3) = T(1:3,1:3);
x(4:6,4:6) = T(1:3,1:3);
x(4:6,1:3) = PSkew * T(1:3,1:3);



end
function x = SO3(angles)% convert ZYX Euler angle to an orientation matrix
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