% Jamison McGinley, Jarrod Puseman
% 4/10/20
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function ddt = odefun2(t,IC)
xi_g = IC(1);
eta_g = IC(2);
theta_g = IC(3);
xi_a = IC(4);
eta_a = IC(5);
theta_a = IC(6);

vg = 2; %m/s
L = 0.5; %m
phi_g = -pi/18; %rad
va = 12; %m/s
wa = pi/25; %rad/s

%Changes in states

dxi_g = vg*cos(theta_g);
deta_g = vg*sin(theta_g);
dtheta_g = (vg/L)*tan(phi_g);
dxi_a = va*cos(theta_a);
deta_a = va*cos(theta_a);
dtheta_a = wa;

%Measurements
dy1 = atan((eta_a - eta_g)/(xi_a - xi_g))-theta_g;
dy2 = sqrt((eta_g - eta_a)^2 + (xi_g - xi_a)^2);
dy3 = atan((eta_g - eta_a)/(xi_g - xi_a))-theta_a;
dy4 = xi_a;
dy5 = eta_a;

ddt = [dxi_g;deta_g;dtheta_g;dxi_a;deta_a;dtheta_a;dy1;dy2;dy3;dy4;dy5];

end