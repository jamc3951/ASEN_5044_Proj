function fig = filterEstimate(mu,P,t,stateNames,plottitle,val)
% filterEstimate plots the estimate of a dynamical system xhat with iteration
% along with 2-sigma bounds extracted from the record of covariance
% matrices P
% Format of call filterEstimate(xhat,P,t,stateNames,units,plot_title,val)
% Returns fig a figure handle to the plot.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ASEN 5044: Statistical Estimation of Dynamic Systems
% Final Project
% Jamison McGinley, Jarrod Puseman
% Dr. Matsuo
% 5/1/2020
% Created:  4/23/2020
% Modified: 4/23/2020
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fig = figure;
lw = 1;
for i  = 1:6% Assumes 6 state Variables
    subplot(2,3,i)
    plot(t,mu(i,:),'b','LineWidth',lw);
    hold on
    grid on
    grid minor
    sigma =  sqrt(reshape(P(i,i,:),1,length(mu)));
    plot(t,val*mu(i,:)+2*sigma,'r--','LineWidth',lw)
    plot(t,val*mu(i,:)-2*sigma,'r--','LineWidth',lw)
    xlabel('Time [s]','interpreter','latex')
    ylabel([stateNames{i}])
    title(stateNames{i})
    if i == 1
        legend(['Estimated ' stateNames{i}],['True ' stateNames{i}],'2-\sigma Error Bounds')
    end
end
suptitle(plottitle)
set(gcf, 'Position', [100, 100, 1100, 730]) %Reposition
end