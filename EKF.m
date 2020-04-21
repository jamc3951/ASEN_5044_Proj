function [x,P,NIS,innovation] = EKF(x0,P0,time,f,A,Omegak,Q,R,h,C,y)
% EKF implements an extended Kalman Filter
% Format of call: 
% Returns: [x,P,NIS] where x is a matrix with the state estimates as a
% function of time and P is a 3D matrix giving the covariance at each time
% and eps is the normalized innovation squared for all k
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ASEN 5044: Statistical Estimation of Dynamic Systems
% Final Project
% Jamison McGinley, Jarrod Puseman
% 5/1/2020
% Created:  4/17/2020
% Modified: 4/17/2020
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
n = length(x0);
len = length(time);
I = eye(n);
p = length(R);

x = zeros(n,len);
P = zeros(n,n,len);
NIS = zeros(1,len-1);
innovation = zeros(p,len-1);
x(:,1) = x0;
P(:,:,1) = P0;
covaradd = Omegak*Q*Omegak';
t = time(1);
dt = time(2)-time(1);
for i = 2:len;
    u = x(:,i-1);
    k1 = dt*f(t,u);
    k2 = dt*f(t+dt/2,u+k1/2);
    k3 = dt*f(t+dt/2,u+k2/2);
    k4 = dt*f(t+dt,u+k3);
    xminus = u+(k1+2*k2+2*k3+k4)/6; %Use RK4 to numerically approx.
    
    F = (I + dt*A(x(:,i-1)));
    Pminus = F*P(:,:,i-1)*F' + covaradd;
       
    H = C(xminus);
    Sk = H*Pminus*H' + R;
    Sk = .5*(Sk+Sk'); %Make positive definite?
    yhat = h(xminus); %Full nonlinear dynamics
    K = Pminus*H'/Sk; % MATLAB will take inverse
    eyk = y(:,i-1) - yhat;
    innovation(:,i-1) = eyk; 
    NIS(i-1) = eyk'/Sk*eyk; %Normalized Innovation Squared (NIS)
    x(:,i) = xminus + K*(eyk);
    P(:,:,i) = (I-K*H)*Pminus;
end
end