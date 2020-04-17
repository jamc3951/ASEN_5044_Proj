%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ASEN 5044: Statistical Estimation of Dynamic Systems
% Final Project
% Jamison McGinley, Jarrod Puseman
% Dr. Matsuo
% 5/1/2020
% Created:  4/10/2020
% Modified: 4/16/2020
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Setup
clear; close all; clc;

% Make plots?
plotbool = [0 1];

dt = 0.1;
vg = 2; %m/s
L = 0.5; %m
phi_g = -pi/18; %rad
va = 12; %m/s
wa = pi/25; %rad/s

n=6;
p=5;
d=4;

xi_g_nom =@(t) (1/(2*tan(pi/18)))*(20*tan(pi/18) + 1 - cos(4*tan(phi_g)*t)); 
eta_g_nom =@(t) (1/(2*tan(pi/18)))*sin(4*tan(-phi_g)*t);
theta_g_nom =@(t) wrapToPi(pi/2 + 4*tan(phi_g)*t);
xi_a_nom =@(t) (1/pi)*(300 - 60*pi - 300*cos(pi/25*t));
eta_a_nom =@(t) -(300/pi)*sin(pi/25*t);
theta_a_nom =@(t) wrapToPi(-pi/2 + pi/25*t);
nom_cond =@(t) [xi_g_nom(t); eta_g_nom(t); theta_g_nom(t); xi_a_nom(t); eta_a_nom(t); theta_a_nom(t)];
perturbation = [0.15;0.15;0.05;0.15;0.15;0.05]; 
%perturbation = [0;1;0; 0;0;0.1]; %Used to compare to TA solution
inishcondish = nom_cond(0);
perturbed_state = inishcondish + perturbation;

meas =@(states) [wrapToPi(atan2((states(5,:) - states(2,:)),(states(4,:) - states(1,:)))-states(3,:)); ...
sqrt((states(2,:) - states(5,:)).^2 + (states(1,:) - states(4,:)).^2); ...
wrapToPi(atan2((states(2,:) - states(5,:)),(states(1,:) - states(4,:)))-states(6,:)); ...
states(4,:); ...
states(5,:)];

load('cooplocalization_finalproj_KFdata.mat');
% measlabels
% Qtrue
% Rtrue
% tvec
% ydata

%% Use ode45 To Predict State
[t_ode, x_ode] = ode45(@(t,y) odefun2(t,y,vg,L,phi_g,va,wa,zeros(n,n)), [0 100], perturbed_state);
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
xnom = nom_cond(time);
ynom = meas(xnom(:,2:end));

for k=2:len-1
    Fk = (eye(6) + dt*A(xnom(:,k-1)));
    Gk = dt*B(time(k-1));
    Omegak = dt*Gamma;
    
    dx(:,k) = Fk*dx(:,k-1) + Gk*du + Omegak*w(:,k);
    dy(:,k) = C(xnom(:,k))*dx(:,k)+ v(:,k);
end
x = xnom+dx;
y = ynom+dy;

%% Plot and compare the two formulations to verify Dynamics
% UGV States
ugvstates = {'\xi_g [m]','\eta_g [m]','\theta_g [rad]'};
uavstates = {'\xi_a [m]','\eta_a [m]','\theta_a [rad]'};

if plotbool(1)
    plotcompare(t_ode,x_ode(:,1:3)',time,x(1:3,:),ugvstates,'Linearized vs. Nonlinear UGV States');
    print('UGV','-dpng')

    plotcompare(t_ode,x_ode(:,4:6)',time,x(4:6,:),uavstates,'Linearized vs. Nonlinear UAV States');
    print('UAV','-dpng')

    plotcompare(t_ode(2:end),measurements,time(2:end),y,measLabels,'Linearized vs. Nonlinear Measurements');
    print('measurements','-dpng')
end

%% Linearized KF - Compute what we can offline
FkLKF = zeros(n,n,len-1); % F(1) is F_0
GkLKF = zeros(n,d,len-1); % G(1) is G_0
OmegakLKF = dt*Gamma; %Constant for all time
HkLKF = zeros(p,n,len-1); %H(1) is H_1
dukLKF = zeros(d,len-1); % No force deviation from nominal

%Fill in the Jacobians
for i = 1:len-1
    FkLKF(:,:,i) = (eye(6) + dt*A(xnom(:,i)));
    GkLKF(:,:,i) = dt*B(time(i));
    HkLKF(:,:,i) = C(xnom(:,i+1)); %Recall indexing different
end

% Monte Carlo
num = 5;
% Make some structure to store all the trials of data
P0 = diag([.25,.25,.1,3,3,.1].^2);
NEES = zeros(num,len);
NIS = zeros(num,len-1);

QEKF = Qtrue;

for i = 1:num
    % Generate x0 
    x0 = mvnrnd(inishcondish,P0)';
    
    % Generate some truth data 
%     % Using ode45
%     [time_ode,xMC2] = ode45(@(t,y) odefun2(t,y,vg,L,phi_g,va,wa,Qtrue), time, x0);
%     xMC2 = xMC2';
%     xMC2([3 6],:) = wrapToPi(xMC2([3 6],:));
    % Using RK4 - Basically the same as ode45 method
    xMC = zeros(6,len);
    xMC(:,1) = x0;
    u = x0;
    t = time(1);
    for j = 2:len
        k1 = dt*odefun2(t,u,vg,L,phi_g,va,wa,zeros(n,n));
        k2 = dt*odefun2(t+dt/2,u+k1/2,vg,L,phi_g,va,wa,zeros(n,n));
        k3 = dt*odefun2(t+dt/2,u+k2/2,vg,L,phi_g,va,wa,zeros(n,n));
        k4 = dt*odefun2(t+dt,u+k3,vg,L,phi_g,va,wa,zeros(n,n));
        u = u+(k1+2*k2+2*k3+k4)/6 + dt*mvnrnd(zeros(1,6),Qtrue)';
        xMC(:,j)= u;
        t = t+dt;
    end
    xMC([3 6],:) = wrapToPi(xMC([3 6],:));
    % Using linearization
%     dx = zeros(n,len);
%     dx(:,1) = x0-inishcondish;
%     w = mvnrnd(zeros(1,6),Qtrue,len)';
%     for k=2:len-1
%         Fk = FkLKF(:,:,k-1);
%         Gk = GkLKF(:,:,k-1);
%         
%         dx(:,k) = Fk*dx(:,k-1) + Gk*du + Omegak*w(:,k);
%     end
%     xMC3 = xnom + dx;
%     %Plot to compare the methods
%     figure
%     titles = {'ode45 with error','RK4','Linearization'};
%     data = {xMC2,xMC,xMC3};
%     for k=1:3
%         subplot(1,3,k)
%         plot(time,data{k}(1,:))
%         title(titles{k})
%         ylabel('\xi_g [m]')
%         xlabel('Time [s]')
%         ylim([9 17])
%     end
    
    
    % Generate measurements 
    measMC = meas(xMC(:,2:end)); %Exact measurements for the states
    noisymeas = measMC + mvnrnd([0;0;0;0;0],Rtrue,len-1)'; %added noise
    % Create dy 
    dyLKF =noisymeas-ynom; %Difference between measurement and nomimal meas
    
    %Linearized KF
    [dxLKF,P,NIS(i,:)] = LKF(zeros(6,1),P0,time,FkLKF,GkLKF,dukLKF,OmegakLKF,QEKF,Rtrue,HkLKF,dyLKF);
    xLKF = xnom + dxLKF;
    
    % Compute NEES
    epsx = xMC-xLKF;
    for j = 1:len
        NEES(i,j) = epsx(:,j)'*(P(:,:,j)\epsx(:,j));
    end
end
NEES_avg = mean(NEES);
NIS_avg = mean(NIS);

% Make the plots requested (a single trial - use xLKF and P, NEES and NIS)
alpha = .05; %significance level
r1 = chi2inv(alpha/2,num*n)./num;
r2 = chi2inv(1-alpha/2,num*n)./num;
r3 = chi2inv(alpha/2,num*p)./num;
r4 = chi2inv(1-alpha/2,num*p)./num;

if plotbool(2)
    %Plot the states
    plotEstimate(xLKF,P,time,xMC,[ugvstates uavstates],'Comparing Truth and Predicted States for a Trial',1);
    print('LKF_est_tuning','-dpng')
    
    %Plot the delta_states
    plotEstimate(dxLKF,P,time,(xMC-xnom),[ugvstates uavstates],'Comparing Delta States for a Trial',0);
    print('LKF_delta_est_tuning','-dpng')
    
    %Plot the errors
    plotEstimate(xMC-xLKF,P,time,zeros(n,len),[ugvstates uavstates],'State Estimation Errors',0)
    print('LKF_error_tuning','-dpng')
    
    figure
    subplot(1,2,1)
    plot(time,NEES_avg,'ko','LineWidth',1)
    hold on
    plot([0,time(end)],[r1 r1],'--r','LineWidth',2)
    plot([0,time(end)],[r2 r2],'--r','LineWidth',2)
    title('NEES Estimation Results')
    ylabel('NEES Statistic, $\bar{\epsilon}_x$','interpreter','latex','Fontsize',14)
    xlabel('Time [s]')
    
    subplot(1,2,2)
    plot(time(2:end),NIS_avg,'ko','LineWidth',1)
    hold on
    plot([time(2),time(end)],[r3 r3],'--r','LineWidth',2)
    plot([time(2),time(end)],[r4 r4],'--r','LineWidth',2)
    title('NIS Estimation Results')
    ylabel('NIS Statistic, $\bar{\epsilon}_y$','interpreter','latex','Fontsize',14)
    xlabel('Time [s]')
    set(gcf, 'Position', [100, 100, 1100, 730]) 
    print('LKF_NEESNIS','-dpng')
end

%% Extended Kalman Filter
% Repeat the above









%% Implement Filters on the Provided Data
% Hopefully easy


