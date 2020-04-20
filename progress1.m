%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ASEN 5044: Statistical Estimation of Dynamic Systems
% Final Project
% Jamison McGinley, Jarrod Puseman
% Dr. Matsuo
% 5/1/2020
% Created:  4/10/2020
% Modified: 4/17/2020
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Setup
clear; close all; clc;

% Make plots?
plotbool = [0 1 1];
rng(100);

dt = 0.1;
vg = 2; %m/s
L = 0.5; %m
phi_g = -pi/18; %rad
va = 12; %m/s
wa = pi/25; %rad/s

n=6;
p=5;
d=4;

f =@(t,x) [vg*cos(x(3));vg*sin(x(3));(vg/L)*tan(phi_g);va*cos(x(6));va*sin(x(6));wa];

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
% [t_ode, x_ode] = ode45(@(t,y) odefun2(t,y,vg,L,phi_g,va,wa), [0 100], perturbed_state);
[t_ode, x_ode] = ode45(f, [0 100], perturbed_state);
measurements = meas(x_ode(2:end,:)');
x_ode(:,[3 6]) = wrapToPi(x_ode(:,[3 6]));

%% Simulate the Linearized System w/ same perturbation
A =@(x) [0 0 -vg*sin(x(3)) 0 0 0; 0 0 vg*cos(x(3)) 0 0 0; ...
     0 0 0 0 0 0; 0 0 0 0 0 -va*sin(x(6)); 0 0 0 0 0 va*cos(x(6)); ...
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

for k=2:len
    Fk = (eye(6) + dt*A(xnom(:,k-1)));
    Gk = dt*B(time(k-1));
    Omegak = dt*Gamma;
    
    dx(:,k) = Fk*dx(:,k-1) + Gk*du + Omegak*w(:,k-1);
    dy(:,k-1) = C(xnom(:,k))*dx(:,k)+ v(:,k-1);
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


%% Monte Carlo - Generate some trials
num = 50;
%Generate some truth data sets
P0 = 100*diag([.5,.5,.1,2,2,.5].^2);
xMC = zeros(n,len,num); %Create same trial runs for each filter
noisymeas = zeros(p,len-1,num);
for i = 1:num
    % Generate x0 
    x0 = mvnrnd(inishcondish,P0)';
    
    % Generate some truth data 
    % Using ode45
%     xMC = zeros(n,len);
%     xMC(:,1) = x0;
%     for j = 2:len
%         [~,statespred] = ode45(@(t,y) odefun2(t,y,vg,L,phi_g,va,wa),[time(j-1) time(j)], xMC(:,j-1));
%         xMC(:,j) = statespred(end,:)' + mvnrnd(zeros(1,6),OmegakLKF*Qtrue*OmegakLKF')';
%     end
%     xMC([3 6],:) = wrapToPi(xMC([3 6],:));

    % Using RK4 - Basically the same as ode45 method, but faster
    xMC(:,1,i) = x0;
    u = x0;
    t = time(1);
    for j = 2:len
        k1 = dt*f(t,u);
        k2 = dt*f(t+dt/2,u+k1/2);
        k3 = dt*f(t+dt/2,u+k2/2);
        k4 = dt*f(t+dt,u+k3);
        u = u+(k1+2*k2+2*k3+k4)/6 + mvnrnd(zeros(1,6),Omegak*Qtrue*Omegak')';
        xMC(:,j,i)= u;
        t = t+dt;
    end
    xMC([3 6],:,i) = wrapToPi(xMC([3 6],:,i));    
    
    % Generate measurements 
    measMC = meas(xMC(:,2:end,i)); %Exact measurements for the states
    noisymeas(:,:,i) = measMC + mvnrnd([0;0;0;0;0],Rtrue,len-1)'; %added noise
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

NEES_LKF = zeros(num,len);
NIS_LKF = zeros(num,len-1);

QLKF = 50000*Qtrue;
RLKF = 80*Rtrue;
for i=1:num
    % Create dy 
    dyLKF =noisymeas(:,:,i)-ynom; %Difference between measurement and nomimal meas
    
    %Linearized KF
    [dxLKF,P,NIS_LKF(i,:)] = LKF(zeros(6,1),P0,time,FkLKF,GkLKF,dukLKF,OmegakLKF,QLKF,RLKF,HkLKF,dyLKF);
    xLKF = xnom + dxLKF;
    
    % Compute NEES
    epsx = xMC(:,:,i)-xLKF;
    for j = 1:len
        NEES_LKF(i,j) = epsx(:,j)'*(P(:,:,j)\epsx(:,j));
    end
end
NEES_avg_LKF = mean(NEES_LKF);
NIS_avg_LKF = mean(NIS_LKF);

% Make the plots requested (a single trial - use xLKF and P, NEES and NIS)
alpha = .05; %significance level
r1 = chi2inv(alpha/2,num*n)./num;
r2 = chi2inv(1-alpha/2,num*n)./num;
r3 = chi2inv(alpha/2,num*p)./num;
r4 = chi2inv(1-alpha/2,num*p)./num;

if plotbool(2)
    %Plot the states
    plotEstimate(xLKF,P,time,xMC(:,:,num),[ugvstates uavstates],'Comparing Truth and Predicted States for a Trial, LKF',1);
    print('LKF_est_tuning','-dpng')
    
    %Plot the delta_states
    dev = xMC(:,:,num)-xnom;
    dev([3 6],:) = wrapToPi(dev([3 6],:));
    plotEstimate(dxLKF,P,time,dev,[ugvstates uavstates],'Comparing Delta States for a Trial, LKF',1);
    print('LKF_delta_est_tuning','-dpng')
    
    %Plot the errors
    err = xMC(:,:,num)-xLKF;
    err([3 6],:) = wrapToPi(err([3 6],:));
    plotEstimate(err,P,time,zeros(n,len),[ugvstates uavstates],'State Estimation Errors, LKF',0);
    print('LKF_error_tuning','-dpng')
    
    figure
    subplot(1,2,1)
    plot(time,NEES_avg_LKF,'ko','LineWidth',1)
    hold on
    plot([0,time(end)],[r1 r1],'--r','LineWidth',2)
    plot([0,time(end)],[r2 r2],'--r','LineWidth',2)
    grid on
    grid minor
    title('NEES Estimation Results, LKF')
    ylabel('NEES Statistic, $\bar{\epsilon}_x$','interpreter','latex','Fontsize',14)
    xlabel('Time [s]')
    
    subplot(1,2,2)
    plot(time(2:end),NIS_avg_LKF,'ko','LineWidth',1)
    hold on
    plot([time(2),time(end)],[r3 r3],'--r','LineWidth',2)
    plot([time(2),time(end)],[r4 r4],'--r','LineWidth',2)
    grid on 
    grid minor
    title('NIS Estimation Results, LKF')
    ylabel('NIS Statistic, $\bar{\epsilon}_y$','interpreter','latex','Fontsize',14)
    xlabel('Time [s]')
    set(gcf, 'Position', [100, 100, 1100, 730]) 
    suptitle('Linearized KF \chi^2 Statistics')
    print('LKF_NEESNIS','-dpng')
end

%% Extended Kalman Filter
NEES_EKF = zeros(num,len);
NIS_EKF = zeros(num,len-1);

REKF = Rtrue;
QEKF = 20000*Qtrue;
QEKF(4:6,:) = 20000*Qtrue(4:6,:);
for i=1:num  
    %Extended KF
    [xEKF,P,NIS_EKF(i,:)] = EKF(inishcondish,P0,time,f,A,Omegak,QEKF,REKF,meas,C,noisymeas(:,:,i));
    xEKF([3 6],:) = wrapToPi(xEKF([3 6],:));
    
    % Compute NEES
    epsx = xMC(:,:,i)-xEKF;
    for j = 1:len
        NEES_EKF(i,j) = epsx(:,j)'*(P(:,:,j)\epsx(:,j));
    end
end
NEES_avg_EKF = mean(NEES_EKF);
NIS_avg_EKF = mean(NIS_EKF);

if plotbool(3)
    %Plot the states
    plotEstimate(xEKF,P,time,xMC(:,:,num),[ugvstates uavstates],'Comparing Truth and Predicted States for a Trial, EKF',1);
    print('EKF_est_tuning','-dpng')
        
    %Plot the errors
    err = xMC(:,:,num)-xEKF;
    err([3 6],:) = wrapToPi(err([3 6],:));
    plotEstimate(err,P,time,zeros(n,len),[ugvstates uavstates],'State Estimation Errors, EKF',0);
    print('EKF_error_tuning','-dpng')
    
    figure
    subplot(1,2,1)
    plot(time,NEES_avg_EKF,'ko','LineWidth',1)
    hold on
    plot([0,time(end)],[r1 r1],'--r','LineWidth',2)
    plot([0,time(end)],[r2 r2],'--r','LineWidth',2)
    grid on
    grid minor
    title('NEES Estimation Results, EKF')
    ylabel('NEES Statistic, $\bar{\epsilon}_x$','interpreter','latex','Fontsize',14)
    xlabel('Time [s]')
    
    subplot(1,2,2)
    plot(time(2:end),NIS_avg_EKF,'ko','LineWidth',1)
    hold on
    plot([time(2),time(end)],[r3 r3],'--r','LineWidth',2)
    plot([time(2),time(end)],[r4 r4],'--r','LineWidth',2)
    grid on
    grid minor
    title('NIS Estimation Results, EKF')
    ylabel('NIS Statistic, $\bar{\epsilon}_y$','interpreter','latex','Fontsize',14)
    xlabel('Time [s]')
    set(gcf, 'Position', [100, 100, 1100, 730]) 
    suptitle('Extended KF \chi^2 Statistics')
    print('EKF_NEESNIS','-dpng')
end


%% Implement Filters on the Provided Data
% Hopefully easy (LOL)


