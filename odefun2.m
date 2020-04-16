%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ASEN 5044: Statistical Estimation of Dynamic Systems
% Final Project
% Jamison McGinley, Jarrod Puseman
% Dr. Matsuo
% 5/1/2020
% Created:  4/10/2020
% Modified: 4/16/2020
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function ddt = odefun2(~,IC,vg,L,phi_g,va,wa,Q)
%xi_g = IC(1);
%eta_g = IC(2);
theta_g = IC(3);
%xi_a = IC(4);
%eta_a = IC(5);
theta_a = IC(6);

w = mvnrnd([0;0;0;0;0;0],Q)'; %Random noise

%Changes in states
dxi_g = vg*cos(theta_g);
deta_g = vg*sin(theta_g);
dtheta_g = (vg/L)*tan(phi_g);
dxi_a = va*cos(theta_a);
deta_a = va*sin(theta_a);
dtheta_a = wa;

ddt = [dxi_g;deta_g;dtheta_g;dxi_a;deta_a;dtheta_a]+w;
end