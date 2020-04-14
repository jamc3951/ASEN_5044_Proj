% Jamison McGinley, Jarrod Puseman
% 4/10/20
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function ddt = odefun2(~,IC)
%xi_g = IC(1);
%eta_g = IC(2);
theta_g = IC(3);
%xi_a = IC(4);
%eta_a = IC(5);
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
deta_a = va*sin(theta_a);
dtheta_a = wa;

ddt = [dxi_g;deta_g;dtheta_g;dxi_a;deta_a;dtheta_a];
end