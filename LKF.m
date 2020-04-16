function [dx,P,eps] = LKF(dx0,P0,time,Fk,Gk,du,Omegak,Q,R,Hk,dy)
% LKF implements a linearized Kalman Filter
% Format of call: 
% Returns: [dx,P,eps] where dx is a matrix with the state estimates as a
% function of time and P is a 3D matrix giving the covariance at each time
% and eps is the normalized innovation squared for all k
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ASEN 5044: Statistical Estimation of Dynamic Systems
% Final Project
% Jamison McGinley, Jarrod Puseman
% 5/1/2020
% Created:  4/16/2020
% Modified: 4/16/2020
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
n = length(dx0);
len = length(time);
I = eye(n);

dx = zeros(n,len);
P = zeros(n,n,len);
eps = zeros(1,len-1);
dx(:,1) = dx0;
P(:,:,1) = P0;
for i = 2:len;
    F = Fk(:,:,i-1);
    dxminus = F*dx(:,i-1) + Gk(:,:,i-1)*du(:,i-1);
    Pminus = F*P(:,:,i-1)*F' + Omegak*Q*Omegak';
       
    H = Hk(:,:,i-1);
    K = Pminus*H'*(H*Pminus*H' + R)^-1;
    dx(:,i) = dxminus + K*(dy(:,i-1) - H*dxminus);
    P(:,:,i) = (I-K*H)*Pminus;
    
    Sk = H*Pminus*H' + R;
    eyk = dy(:,i-1)-H*dxminus;
    eps(i-1) = eyk'*(Sk\eyk);
end
end