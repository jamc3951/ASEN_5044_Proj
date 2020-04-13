clear; close all; clc;

% Use ode45 to get out expected results from the Nonlinear System

dt = 0.1;
vg = 2; %m/s
L = 0.5; %m
phi_g = -pi/18; %rad
va = 12; %m/s
wa = pi/25; %rad/s

xi_g_nom =@(t) (1/(2*tan(pi/18)))*(20*tan(pi/18) + 1 - cos(4*tan(-pi/18)*t)); 
eta_g_nom =@(t) (1/(2*tan(pi/18)))*sin(4*tan(pi/18)*t);
theta_g_nom =@(t) pi/2 - 4*tan(pi/18)*t;
xi_a_nom =@(t) (1/pi)*(300 - 60*pi - 300*cos(pi/25*t));
eta_a_nom =@(t) -(300/pi)*sin(pi/25*t);
theta_a_nom =@(t) -pi/2 + pi/25*t;

nom_cond =@(t) [xi_g_nom(t); eta_g_nom(t); theta_g_nom(t); xi_a_nom(t); eta_a_nom(t); theta_a_nom(t)];
perturbation = [0.15,0.15,0.05,0.15,0.15,0.05]; 
%perturbation = zeros(1,6);
inishcondish = [10,0,pi/2,-60,0,-pi/2];
perturbed_state = inishcondish + perturbation;


[t_ode, x_ode] = ode45('odefun2', [0 100], perturbed_state);
xi_g = x_ode(2:end,1);
eta_g = x_ode(2:end,2);
theta_g = x_ode(2:end,3);
xi_a = x_ode(2:end,4);
eta_a = x_ode(2:end,5);
theta_a = x_ode(2:end,6);

meas =@(states) [atan2((states(:,5) - states(:,2)),(states(:,4) - states(:,1)))-states(:,3), ...
sqrt((states(:,2) - states(:,5)).^2 + (states(:,1) - states(:,4)).^2), ...
atan2((states(:,2) - states(:,5)),(states(:,1) - states(:,4)))-states(:,6), ...
states(:,4), ...
states(:,5)];
measurements = meas(x_ode(2:end,:));
% Simulate the Linearized System w/ same perturbation


A =@(t) [0 0 -2*sin(pi/2 + 4*tan(pi/18)*t) 0 0 0; 0 0 2*cos(pi/2 + 4*tan(pi/18)*t) 0 0 0; ...
     0 0 0 0 0 0; 0 0 0 0 0 -12*sin(-pi/2 + pi/25*t); 0 0 0 0 0 12*cos(-pi/2 + pi/25*t); ...
     0 0 0 0 0 0];

time = 0:dt:100;
len = length(time);
dx = zeros(6,len);
dy = zeros(5,len-1);
du = [0;0;0;0];
dx(:,1)= perturbation;

for k=2:len-1
    
    Fk = (eye(6) + dt*A(time(k-1)));
    Gk = zeros(6,4);
    dx(:,k) = Fk*dx(:,k-1) + Gk*du;
    
    nomstate = nom_cond(time(k));
    xi_g = nomstate(1);
    eta_g = nomstate(2);
    xi_a = nomstate(4);
    eta_a = nomstate(5);
    
    
    C = [(eta_a - eta_g)/(((eta_g - eta_a)^2 + (xi_g - xi_a)^2)), (-xi_a + xi_g)/(((eta_g - eta_a)^2 + (xi_g - xi_a)^2)), -1, (-eta_a + eta_g)/(((eta_g - eta_a)^2 + (xi_g - xi_a)^2)), (xi_a - xi_g)/(((eta_g - eta_a)^2 + (xi_g - xi_a)^2)), 0; ...
    (xi_g - xi_a)/sqrt((eta_g - eta_a)^2 + (xi_g - xi_a)^2), (eta_g - eta_a)/sqrt((eta_g - eta_a)^2 + (xi_g - xi_a)^2), 0, -(xi_g - xi_a)/sqrt((eta_g - eta_a)^2 + (xi_g - xi_a)^2), -(eta_g - eta_a)/sqrt((eta_g - eta_a)^2 + (xi_g - xi_a)^2), 0; ...
    (eta_a - eta_g)/(((eta_g - eta_a)^2 + (xi_g - xi_a)^2)), (-xi_a + xi_g)/(((eta_g - eta_a)^2 + (xi_g - xi_a)^2)), 0, (-eta_a + eta_g)/(((eta_g - eta_a)^2 + (xi_g - xi_a)^2)), (-xi_a + xi_g)/(((eta_g - eta_a)^2 + (xi_g - xi_a)^2)), -1; ...
    0 0 0 1 0 0; 0 0 0 0 1 0];
    
    dy(:,k) = C*dx(:,k);
    
    
end
x = nom_cond(time)+dx;
y = meas(nom_cond(time(2:end))')+dy';
%Plot and compare the two formulations

figure()
hold on;
grid on;
plot(t_ode,x_ode(:,[1,2,4,5]));
plot(time,x([1,2,4,5],:));

figure()
hold on;
grid on;
plot(t_ode(2:end),measurements(:,:));
plot(time(2:end),y(:,:));




















