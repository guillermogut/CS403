function u = hw7controller(z, param, x_des, dx_des, ddx_des) 
   % ******** Implement your controller ********
   
   key_pt = keypoints_pend(z,param);
   rB = key_pt(:,2);
   vB = velocity_rB(z,param);
   
   errPos =  x_des - rB;
   errVel = dx_des - vB;
   Jb = Jacobian_rB(z,param);
   Jbdot = Jdot_rB(z,param);
   C = coriolis_pend(z,param);
   M = a_pend(z,param);
   G = grav_pend(z,param);
   lambdaq = inv(Jb*pinv(M)*Jb.');
   mu = (lambdaq * Jb * inv(M)*C)- lambdaq*Jbdot*z(3:4);
   rho = lambdaq*Jb*inv(M)* G;
 
   K = [50 0;
        0 50];
   
   D = [5 0;
        0 5];
  
   %u = Jb.'*(lambdaq*(ddx_des + K*errPos + D*errVel) +mu +rho); 
   %u = Jb.'*(lambdaq*(ddx_des + K*errPos + D*errVel) +mu);% the xy actual seems to lag behind the xy desired but y seems to lag much more
   u = Jb.'*(lambdaq*(ddx_des + K*errPos + D*errVel) +rho);%mu 
   %u = Jb.'*(lambdaq*(K*errPos + D*errVel) +mu +rho);%
   
   
end