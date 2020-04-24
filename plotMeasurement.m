function fig = plotMeasurement(y,Sk,t,truth,stateNames,plottitle,val)
% plotMeasurement plots the measurements of a dynamical system xhat with iteration
% along with 2-sigma bounds extracted from the record of covariance
% matrices Sk
% Format of call plotMeasurement(y,Sk,t,truth,stateNames,units,plot_title,val)
% Returns fig a figure handle to the plot.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ASEN 5044: Statistical Estimation of Dynamic Systems
% Final Project
% Jamison McGinley, Jarrod Puseman
% Dr. Matsuo
% 5/1/2020
% Created:  4/16/2020
% Modified: 4/17/2020
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fig = figure;
lw = 1;
for i  = 1:5% Assumes 6 state Variables
    subplot(5,1,i)
    plot(t,y(i,:),'b','LineWidth',lw);
    hold on
    plot(t,truth(i,:),'k','Linewidth',lw);
    grid on
    grid minor
    sigma =  sqrt(reshape(Sk(i,i,:),1,length(y)));
    plot(t,val*y(i,:)+2*sigma,'r--','LineWidth',lw)
    plot(t,val*y(i,:)-2*sigma,'r--','LineWidth',lw)
    ylabel([stateNames{i}])
    %title(stateNames{i})
    if i == 1
        legend('Estimate ','Truth','2-\sigma Error Bounds')
    end
    if i == 5
        xlabel('Time [s]','interpreter','latex')
    end
end
suptitle(plottitle)
set(gcf, 'Position', [100, 100, 1100, 730]) %Reposition
end