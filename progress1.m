% Jamison McGinley, Jarrod Puseman
% 4/14/20
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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
theta_g_nom =@(t) pi/2 + 4*tan(phi_g)*t;
xi_a_nom =@(t) (1/pi)*(300 - 60*pi - 300*cos(pi/25*t));
eta_a_nom =@(t) -(300/pi)*sin(pi/25*t);
theta_a_nom =@(t) -pi/2 + pi/25*t;
nom_cond =@(t) [xi_g_nom(t); eta_g_nom(t); theta_g_nom(t); xi_a_nom(t); eta_a_nom(t); theta_a_nom(t)];
perturbation = [0.15,0.15,0.05,0.15,0.15,0.05]; 
inishcondish = [10,0,pi/2,-60,0,-pi/2];
perturbed_state = inishcondish + perturbation;

meas =@(states) [atan2((states(:,5) - states(:,2)),(states(:,4) - states(:,1)))-states(:,3), ...
sqrt((states(:,2) - states(:,5)).^2 + (states(:,1) - states(:,4)).^2), ...
atan2((states(:,2) - states(:,5)),(states(:,1) - states(:,4)))-states(:,6), ...
states(:,4), ...
states(:,5)];

%% Use ode45 To Predict State
[t_ode, x_ode] = ode45('odefun2', [0 100], perturbed_state);
measurements = meas(x_ode(2:end,:));

%% Simulate the Linearized System w/ same perturbation
A =@(t) [0 0 -2*sin(theta_g_nom(t)) 0 0 0; 0 0 2*cos(theta_g_nom(t)) 0 0 0; ...
     0 0 0 0 0 0; 0 0 0 0 0 -12*sin(theta_a_nom(t)); 0 0 0 0 0 12*cos(theta_a_nom(t)); ...
     0 0 0 0 0 0];
 
B =@(t) [cos(theta_g_nom(t)), 0, 0 ,0;...
    sin(theta_g_nom(t)), 0 , 0, 0;...
    2*tan(phi_g), 4*sec(phi_g)^2, 0, 0;...
    0,0, cos(theta_a_nom(t)),0;...
    0,0, sin(theta_a_nom(t)),0;...
    0,0,0,1];

C =@(t) [(eta_a_nom(t) - eta_g_nom(t))/(((eta_g_nom(t) - eta_a_nom(t))^2 + (xi_g_nom(t) - xi_a_nom(t))^2)), (-xi_a_nom(t) + xi_g_nom(t))/(((eta_g_nom(t) - eta_a_nom(t))^2 + (xi_g_nom(t) - xi_a_nom(t))^2)), -1, (-eta_a_nom(t) + eta_g_nom(t))/(((eta_g_nom(t) - eta_a_nom(t))^2 + (xi_g_nom(t) - xi_a_nom(t))^2)), (xi_a_nom(t) - xi_g_nom(t))/(((eta_g_nom(t) - eta_a_nom(t))^2 + (xi_g_nom(t) - xi_a_nom(t))^2)), 0; ...
    (xi_g_nom(t) - xi_a_nom(t))/sqrt((eta_g_nom(t) - eta_a_nom(t))^2 + (xi_g_nom(t) - xi_a_nom(t))^2), (eta_g_nom(t) - eta_a_nom(t))/sqrt((eta_g_nom(t) - eta_a_nom(t))^2 + (xi_g_nom(t) - xi_a_nom(t))^2), 0, -(xi_g_nom(t) - xi_a_nom(t))/sqrt((eta_g_nom(t) - eta_a_nom(t))^2 + (xi_g_nom(t) - xi_a_nom(t))^2), -(eta_g_nom(t) - eta_a_nom(t))/sqrt((eta_g_nom(t) - eta_a_nom(t))^2 + (xi_g_nom(t) - xi_a_nom(t))^2), 0; ...
    (eta_a_nom(t) - eta_g_nom(t))/(((eta_g_nom(t) - eta_a_nom(t))^2 + (xi_g_nom(t) - xi_a_nom(t))^2)), (-xi_a_nom(t) + xi_g_nom(t))/(((eta_g_nom(t) - eta_a_nom(t))^2 + (xi_g_nom(t) - xi_a_nom(t))^2)), 0, (-eta_a_nom(t) + eta_g_nom(t))/(((eta_g_nom(t) - eta_a_nom(t))^2 + (xi_g_nom(t) - xi_a_nom(t))^2)), (-xi_a_nom(t) + xi_g_nom(t))/(((eta_g_nom(t) - eta_a_nom(t))^2 + (xi_g_nom(t) - xi_a_nom(t))^2)), -1; ...
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

for k=2:len-1
    Fk = (eye(6) + dt*A(time(k-1)));
    Gk = dt*B(time(k-1));
    Omegak = dt*Gamma;
    
    dx(:,k) = Fk*dx(:,k-1) + Gk*du + Omegak*w(:,k);
    dy(:,k) = C(time(k))*dx(:,k)+ v(:,k);
end
x = nom_cond(time)+dx;
y = meas(nom_cond(time(2:end))')+dy';

%% Plot and compare the two formulations
% UGV States
figure()
subplot(3,1,1)
hold on;
grid on;
plot(t_ode,x_ode(:,1),'.','Markersize',20);
plot(time,x(1,:),'linewidth',2);
%xlabel('Time [s]','Fontsize',14);
ylabel('\xi_g [m]','Fontsize',14);
legend('ode45', 'Linearized')

subplot(3,1,2)
hold on;
grid on;
plot(t_ode,x_ode(:,2),'.','Markersize',20);
plot(time,x(2,:),'linewidth',2);
%xlabel('Time [s]');
ylabel('\eta_g [m]','Fontsize',14);

subplot(3,1,3)
hold on;
grid on;
plot(t_ode,x_ode(:,3),'.','Markersize',20);
plot(time,x(3,:),'linewidth',2);
xlabel('Time [s]','Fontsize',14);
ylabel('\theta_g [rad]','Fontsize',14);

suptitle('Linearized vs. Nonlinear UGV States')
set(gcf, 'Position', [100, 100, 1100, 730])
print('UGV','-dpng')

% UAV States
figure()
subplot(3,1,1)
hold on;
grid on;
plot(t_ode,x_ode(:,4),'.','Markersize',20);
plot(time,x(4,:),'linewidth',2);
%xlabel('Time [s]');
ylabel('\xi_A [m]','Fontsize',14);
legend('ode45', 'Linearized')

subplot(3,1,2)
hold on;
grid on;
plot(t_ode,x_ode(:,5),'.','Markersize',20);
plot(time,x(5,:),'linewidth',2);
%xlabel('Time [s]');
ylabel('\eta_A [m]','Fontsize',14);


subplot(3,1,3)
hold on;
grid on;
plot(t_ode,x_ode(:,6),'.','Markersize',20);
plot(time,x(6,:),'linewidth',2);
xlabel('Time [s]','Fontsize',14);
ylabel('\theta_A [rad]','Fontsize',14);

suptitle('Linearized vs. Nonlinear UAV States')
set(gcf, 'Position', [100, 100, 1100, 730])
print('UAV','-dpng')


%Measurements
figure()
subplot(5,1,1)
hold on;
grid on;
plot(t_ode(2:end),measurements(:,1),'.','Markersize',20);
plot(time(2:end),y(:,1),'linewidth',2);
%xlabel('Time [s]');
ylabel('y_1 [rad]');
legend('ode45', 'Linearized')

subplot(5,1,2)
hold on;
grid on;
plot(t_ode(2:end),measurements(:,2),'.','Markersize',20);
plot(time(2:end),y(:,2),'linewidth',2);
%xlabel('Time [s]');
ylabel('y_2 [m]');

subplot(5,1,3)
hold on;
grid on;
plot(t_ode(2:end),measurements(:,3),'.','Markersize',20);
plot(time(2:end),y(:,3),'linewidth',2);
%xlabel('Time [s]');
ylabel('y_3 [rad]');

subplot(5,1,4)
hold on;
grid on;
plot(t_ode(2:end),measurements(:,4),'.','Markersize',20);
plot(time(2:end),y(:,4),'linewidth',2);
%xlabel('Time [s]');
ylabel('\xi_A [m]');

subplot(5,1,5)
hold on;
grid on;
plot(t_ode(2:end),measurements(:,5),'.','Markersize',20);
plot(time(2:end),y(:,5),'linewidth',2);
xlabel('Time [s]');
ylabel('\eta_A [m]');

suptitle('Linearized vs. Nonlinear Measurements')
set(gcf, 'Position', [100, 100, 1100, 730])
print('measurements','-dpng')