clear; close; clc;

% Use ode45 to get out expected results from the Nonlinear System

dt = 0.1;
vg = 2; %m/s
L = 0.5; %m
phi_g = -pi/18; %rad
va = 12; %m/s
wa = pi/25; %rad/s

xi_g_nom = (1/(2*tan(pi/18)))*(20*tan(pi/18) - 1 + cos(4*tan(pi/18)*dt));
eta_g_nom = (1/(2*tan(pi/18)))*sin(4*tan(pi/18)*dt);
theta_g_nom = pi/2 + 4*tan(pi/18)*dt;
xi_a_nom = (1/pi)*(300 - 60*pi - 300*cos(pi/25*dt));
eta_a_nom = -(300/pi)*sin(pi/25*dt);
theta_a_nom = -pi/2 + pi/25*dt;

nom_cond = [xi_g_nom, eta_g_nom, theta_g_nom, xi_a_nom, eta_a_nom, theta_a_nom];
perturbation = [0.15,0.15,0.5,0.15,0.15,0.5]; 

perturbed_state = nom_cond + perturbation;

[t_ode, x_ode] = ode45('odefun2', [0 100], [perturbed_state zeros(1,5)]);

states = x_ode(:,1:6);
measurements = x_ode(:,6:11);

% Simulate the Linearized System w/ same perturbation

t = 0.1;
A = [0 0 -2*sin(pi/2 + 4*tan(pi/18)*t) 0 0 0; 0 0 2*cos(pi/2 + 4*tan(pi/18)*t) 0 0 0; ...
     0 0 0 0 0 0; 0 0 0 0 0 -12*sin(-pi/2 + pi/25*t); 0 0 0 0 0 12*cos(-pi/2 + pi/25*t); ...
     0 0 0 0 0 0];

x = zeros(6,1000);
y = zeros(5,1000);
x(:,1) = perturbed_state;

for k=2:200
    x(:,k) = (eye(6) + t*A)*x(:,k-1);
    
    xi_g = x(1,k);
    eta_g = x(2,k);
    xi_a = x(4,k);
    eta_a = x(5,k);
    
    
    C = [(eta_a - eta_g)/(((eta_g - eta_a)^2 + (xi_g - xi_a)^2)), (-xi_a + xi_g)/(((eta_g - eta_a)^2 + (xi_g - xi_a)^2)), -1, (-eta_a + eta_g)/(((eta_g - eta_a)^2 + (xi_g - xi_a)^2)), (xi_a - xi_g)/(((eta_g - eta_a)^2 + (xi_g - xi_a)^2)), 0; ...
    (xi_g - xi_a)/sqrt((eta_g - eta_a)^2 + (xi_g - xi_a)^2), (eta_g - eta_a)/sqrt((eta_g - eta_a)^2 + (xi_g - xi_a)^2), 0, -(xi_g - xi_a)/sqrt((eta_g - eta_a)^2 + (xi_g - xi_a)^2), -(eta_g - eta_a)/sqrt((eta_g - eta_a)^2 + (xi_g - xi_a)^2), 0; ...
    (eta_a - eta_g)/(((eta_g - eta_a)^2 + (xi_g - xi_a)^2)), (-xi_a + xi_g)/(((eta_g - eta_a)^2 + (xi_g - xi_a)^2)), -1, (-eta_a + eta_g)/(((eta_g - eta_a)^2 + (xi_g - xi_a)^2)), (xi_a - xi_g)/(((eta_g - eta_a)^2 + (xi_g - xi_a)^2)), -1; ...
    0 0 0 1 0 0; 0 0 0 0 1 0];
    
    y(:,k) = C*x(:,k);
    
    
end























