%%
addpath('..\matlab_utils')
%syms th1 th2 th3 F1 F2 t1 t2 l1 l2 l3
%

% Question 1 stuff
%eu1 = [0.3, 0.2, 0.5];
eu1 = [0.3, 0.2, 0.5];
eu2 = [0.7, pi, pi/2];
eu3 = [pi/3, 0, 0];



orientation1 = eulerToOrientation(eu1);
orientation2 = eulerToOrientation(eu2);
orientation3 = eulerToOrientation(eu3);


% Question 2 stuff
% find the transformation matrix of frame{3} to frame{0} 0T3 

%transform matrix from frame{0} to frame{3}
% T_03 = T_01*T_12*T_23

P_01 = [0.4, 0.8, 1.2];
R_01 = eulerToOrientation(eu1);
T_01 = SE(R_01,transpose(P_01));

P_12 = [-0.4, 0.5, 1.0];
R_12 = eulerToOrientation(eu2);
T_12 = SE(R_12,transpose(P_12));

P_23 = [0.5, -0.8, 1.2];
R_23 = eulerToOrientation(eu3);
T_23 = SE(R_23,transpose(P_23));

T_02 = T_01 * T_12;
R_02 = T_02(1:3,1:3);
P_02 = T_02(1:3,4);


T_03 = T_01 * T_12 * T_23;
R_03 = T_03(1:3,1:3);
P_03 = T_03(1:3,4);

%Just in case I misunderstood question 2, here is the transformation matrix
%from frame 3 to frame 0

%T_30 = inverse( T_03 ) = inverse(T_01) * inverse(T_12) * inverse(T_23)

T_30 = SEInverse(T_01)*SEInverse(T_12)*SEInverse(T_23);

% Question 3 stuff


% draw figure
figure
hold on
grid on


R = eul2rotm([.0, 0, 0] );
p = [0, 0, 0];

t = 10;
for i=1:30
    
    t = t+0.03;
    drawCoordinate3D(R,p);
    drawCoordinate3D(R_01,P_01);
    
    drawCoordinate3D(R_02,P_02);
    
    drawCoordinate3D(R_03,P_03);
    hold off
  
      pause(.01);
    %xyz vector for the EE
      test = [(0.1*sin(pi*t)+0.05); (0.3*cos(pi*t)+0.08); (sin(pi*t)+0.5);1];
translationVectorForEE = T_03*test;

translationVectorForEE = translationVectorForEE(1:3);



drawCoordinate3D(R_03,translationVectorForEE);
hold off
grid on
xlabel('$x$','interpreter','latex','fontsize',20)
ylabel('$y$','interpreter','latex','fontsize',20)
zlabel('$z$','interpreter','latex','fontsize',20)

axis equal
view(30,30)

end


function x = eulerToOrientation(angles)% convert ZYX Euler angle to an orientation matrix

thZ = angles(1);
thY = angles(2);
thX = angles(3);

Z = [cos(thZ)   -sin(thZ)   0;sin(thZ)     cos(thZ)   0;0           0           1];
Y = [cos(thY)  0  sin(thY);0        1      0; -sin(thY)  0  cos(thY)];
X = [1      0      0;0     cos(thX) -sin(thX);0     sin(thX) cos(thX)];    

x = Z*Y*X;
end    

function se = SE(R,P)%get SE(3) group with rotation R and translation V

leftNoZerosOnes = [R P];

bottom = [0 0 0 1];

se = [leftNoZerosOnes;bottom];

end

function product = mult2SE3(se1,se2)%not sure why we need this one
product = se1*se2;
end

function inverse = SEInverse(m)%inverse transformation matrix 

R = m(1:3,1:3);
P = m(1:3,4);

Rt = transpose(R);
negRt = -1*Rt;

inverse = [Rt negRt*P;0 0 0 1];   
end














