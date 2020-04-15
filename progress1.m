%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ASEN 5044: Statistical Estimation of Dynamic Systems
% Final Project
% Jamison McGinley, Jarrod Puseman
% Dr. Matsuo
% 5/1/2020
% Created:  4/10/2020
% Modified: 4/14/2020
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Setup
clear; close all; clc;

dt = 0.1;
vg = 2; %m/s
L = 0.5; %m
phi_g = -pi/18; %rad
va = 12; %m/s
wa = pi/25; %rad/s

n=6;
p=5;

xi_g_nom =@(t) (1/(2*tan(pi/18)))*(20*tan(pi/18) + 1 - cos(4*tan(phi_g)*t)); 
eta_g_nom =@(t) (1/(2*tan(pi/18)))*sin(4*tan(-phi_g)*t);
theta_g_nom =@(t) wrapToPi(pi/2 + 4*tan(phi_g)*t);
xi_a_nom =@(t) (1/pi)*(300 - 60*pi - 300*cos(pi/25*t));
eta_a_nom =@(t) -(300/pi)*sin(pi/25*t);
theta_a_nom =@(t) wrapToPi(-pi/2 + pi/25*t);
nom_cond =@(t) [xi_g_nom(t); eta_g_nom(t); theta_g_nom(t); xi_a_nom(t); eta_a_nom(t); theta_a_nom(t)];
perturbation = [0.15;0.15;0.05;0.15;0.15;0.05]; 
%perturbation = [0;1;0; 0;0;0.1]; %Used to compare to TA solution
inishcondish = [10;0;pi/2;-60;0;-pi/2];
perturbed_state = inishcondish + perturbation;

meas =@(states) [wrapToPi(atan2((states(5,:) - states(2,:)),(states(4,:) - states(1,:)))-states(3,:)); ...
sqrt((states(2,:) - states(5,:)).^2 + (states(1,:) - states(4,:)).^2); ...
wrapToPi(atan2((states(2,:) - states(5,:)),(states(1,:) - states(4,:)))-states(6,:)); ...
states(4,:); ...
states(5,:)];

%% Use ode45 To Predict State
[t_ode, x_ode] = ode45('odefun2', [0 100], perturbed_state);
measurements = meas(x_ode(2:end,:)');
x_ode(:,[3 6]) = wrapToPi(x_ode(:,[3 6]));

%% Simulate the Linearized System w/ same perturbation
A =@(x) [0 0 -2*sin(x(3)) 0 0 0; 0 0 2*cos(x(3)) 0 0 0; ...
     0 0 0 0 0 0; 0 0 0 0 0 -12*sin(x(6)); 0 0 0 0 0 12*cos(x(6)); ...
     0 0 0 0 0 0];
 
B =@(t) [cos(theta_g_nom(t)), 0, 0 ,0;...
    sin(theta_g_nom(t)), 0 , 0, 0;...
    2*tan(phi_g), 4*sec(phi_g)^2, 0, 0;...
    0,0, cos(theta_a_nom(t)),0;...
    0,0, sin(theta_a_nom(t)),0;...
    0,0,0,1];

C =@(x) [(x(5) - x(2))/(((x(2) - x(5))^2 + (x(1) - x(4))^2)), (-x(4) + x(1))/(((x(2) - x(5))^2 + (x(1) - x(4))^2)), -1, (-x(5) + x(2))/(((x(2) - x(5))^2 + (x(1) - x(4))^2)), (x(4) - x(1))/(((x(2) - x(5))^2 + (x(1) - x(4))^2)), 0; ...
    (x(1) - x(4))/sqrt((x(2) - x(5))^2 + (x(1) - x(4))^2), (x(2) - x(5))/sqrt((x(2) - x(5))^2 + (x(1) - x(4))^2), 0, -(x(1) - x(4))/sqrt((x(2) - x(5))^2 + (x(1) - x(4))^2), -(x(2) - x(5))/sqrt((x(2) - x(5))^2 + (x(1) - x(4))^2), 0; ...
    (x(5) - x(2))/(((x(2) - x(5))^2 + (x(1) - x(4))^2)), (-x(4) + x(1))/(((x(2) - x(5))^2 + (x(1) - x(4))^2)), 0, (-x(5) + x(2))/(((x(2) - x(5))^2 + (x(1) - x(4))^2)), (-x(1) + x(4))/(((x(2) - x(5))^2 + (x(1) - x(4))^2)), -1; ...
    0 0 0 1 0 0; 0 0 0 0 1 0];

Gamma = eye(6);

time = 0:dt:100;
len = length(time);
dx = zeros(n,len);
dy = zeros(p,len-1);
du = [0;0;0;0];
w = zeros(n,len-1); %No noise
v = zeros(p,len-1); %No noise
dx(:,1)= perturbation;
xNom = nom_cond(time);

for k=2:len-1
    Fk = (eye(6) + dt*A(xNom(:,k-1)));
    Gk = dt*B(time(k-1));
    Omegak = dt*Gamma;
    
    dx(:,k) = Fk*dx(:,k-1) + Gk*du + Omegak*w(:,k);
    dy(:,k) = C(xNom(:,k))*dx(:,k)+ v(:,k);
end
x = xNom+dx;
y = meas(xNom(:,2:end))+dy;

%% Plot and compare the two formulations
% UGV States
ugvstates = {'\xi_g [m]','\eta_g [m]','\theta_g [rad]'};
uavstates = {'\xi_a [m]','\eta_a [m]','\theta_a [rad]'};
measurement_labels = {'\gamma_{ag} [rad]','\rho_{ga} [m]','\gamma_{ga} [rad]','\xi_a [m]','\eta_a [m]'};

plotcompare(t_ode,x_ode(:,1:3)',time,x(1:3,:),ugvstates,'Linearized vs. Nonlinear UGV States');
print('UGV','-dpng')

plotcompare(t_ode,x_ode(:,4:6)',time,x(4:6,:),uavstates,'Linearized vs. Nonlinear UAV States');
print('UAV','-dpng')

plotcompare(t_ode(2:end),measurements,time(2:end),y,measurement_labels,'Linearized vs. Nonlinear Measurements');
print('measurements','-dpng')

%% Implement Filters and Test Them

