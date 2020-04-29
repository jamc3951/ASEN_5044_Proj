function [dx,P,innovation,Sk] = LKF(dx0,P0,time,Fk,Gk,du,Omegak,Q,R,Hk,dy)
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
% Modified: 4/25/2020
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
n = length(dx0);
len = length(time);
I = eye(n);
p = length(R);

dx = zeros(n,len);
P = zeros(n,n,len);
innovation = zeros(p,len-1);
Sk = zeros(p,p,len-1);
dx(:,1) = dx0;
P(:,:,1) = P0;
covaradd = Omegak*Q*Omegak';
for i = 2:len;
    F = Fk(:,:,i-1); %Starts at F_0
    dxminus = F*dx(:,i-1) + Gk(:,:,i-1)*du(:,i-1);
    dxminus([3 6]) = wrapToPi(dxminus([3 6]));
    Pminus = F*P(:,:,i-1)*F' + covaradd;
       
    H = Hk(:,:,i-1); %Starts at H_1
    Ski = H*Pminus*(H') + R;
    Ski = .5*(Ski+Ski'); %Make positive definite
    Sk(:,:,i-1) = Ski;
    yhat = H*dxminus;
    yhat([1 3]) = wrapToPi(yhat([1 3]));
    K = (Pminus*(H'))/Ski; % MATLAB will take inverse
    eyk = dy(:,i-1) - yhat;
    %eyk([1 3]) = wrapToPi(eyk([1 3])); %no longer general
    innovation(:,i-1) = eyk;
    dxplus = dxminus + K*(eyk);
    dxplus([3 6]) = wrapToPi(dxplus([3 6]));
    dx(:,i) = dxplus;
    P(:,:,i) = (I-K*H)*Pminus;
end
end