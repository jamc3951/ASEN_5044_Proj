function [dx,P,epsky,innovation] = LKF(dx0,P0,time,Fk,Gk,du,Omegak,Q,R,Hk,dy)
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
p = length(R);

dx = zeros(n,len);
P = zeros(n,n,len);
epsky = zeros(1,len-1);
innovation = zeros(p,len-1);
dx(:,1) = dx0;
P(:,:,1) = P0;
covaradd = Omegak*Q*Omegak';
for i = 2:len;
    F = Fk(:,:,i-1);
    dxminus = F*dx(:,i-1) + Gk(:,:,i-1)*du(:,i-1);
%    dxminus([3 6]) = wrapToPi(dxminus([3 6]));
    Pminus = F*P(:,:,i-1)*F' + covaradd;
       
    H = Hk(:,:,i-1);
    Sk = H*Pminus*H' + R;
    Sk = .5*(Sk+Sk'); %Make positive definite?
    yhat = H*dxminus;
    yhat([1 3]) = wrapToPi(yhat([1 3]));
    K = Pminus*H'/Sk; % MATLAB will take inverse
    eyk = dy(:,i-1) - yhat;
    eyk([1 3]) = wrapToPi(dy([1 3],i-1)-yhat([1 3])); %no longer general
    innovation(:,i-1) = eyk;
    epsky(i-1) = eyk'/Sk*eyk; %Normalized Innovation Squared (NIS)
    dx(:,i) = dxminus + K*(eyk);
    P(:,:,i) = (I-K*H)*Pminus;
end
end