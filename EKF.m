function [x,P,innovation,Sk] = EKF(x0,P0,time,f,A,Omegak,Q,R,h,C,y)
% EKF implements an extended Kalman Filter
% Format of call: 
% Returns: [x,P,innovation,Sk] where x is a matrix with the state estimates as a
% function of time and P is a 3D matrix giving the covariance at each time
% innovation is the innovations in time and their error covariances Sk
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
innovation = zeros(p,len-1);
Sk = zeros(p,p,len-1);
x(:,1) = x0;
P(:,:,1) = P0;
covaradd = Omegak*Q*Omegak';
t = time(1);
dt = time(2)-time(1);
subnum = 30;
subt = dt/subnum;
for i = 2:len;
    u = x(:,i-1); %Not necessary?
    for j = 1:subnum
        k1 = subt*f(t,u);
        k2 = subt*f(t+subt/2,u+k1/2);
        k3 = subt*f(t+subt/2,u+k2/2);
        k4 = subt*f(t+subt,u+k3);
        u = u+(k1+2*k2+2*k3+k4)/6;
        t = t+subt;
    end
    xminus = u; %Use RK4 to numerically approx.

    F = (I + dt*A(x(:,i-1)));
    Pminus = F*P(:,:,i-1)*F' + covaradd;
       
    H = C(xminus);
    Ski = H*Pminus*H' + R;
    Ski = .5*(Ski+Ski'); %Make positive definite?
    Sk(:,:,i-1) = Ski;
    yhat = h(xminus); %Full nonlinear dynamics
    K = Pminus*H'/Ski; % MATLAB will take inverse
    eyk = y(:,i-1) - yhat;
    eyk([1 3]) = wrapToPi(y([1 3],i-1)-yhat([1 3])); %No longer general
    innovation(:,i-1) = eyk; 
    x(:,i) = xminus + K*(eyk);
    P(:,:,i) = (I-K*H)*Pminus;
end
end