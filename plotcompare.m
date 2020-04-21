function fig = plotcompare(tnl,xnl,tl,xl,labels,title)
% Plot Compare creates a figure to compare the propagations of some states
% between a nonlinear solver and linearized solution.
% Format of call: plotcompare(tnl,xnl,tl,xl,lables,title)
% Returns fig a figure handle

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ASEN 5044: Statistical Estimation of Dynamic Systems
% Final Project
% Jamison McGinley, Jarrod Puseman
% Dr. Matsuo
% 5/1/2020
% Created:  4/10/2020
% Modified: 4/14/2020
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

iter = length(labels);
fig = figure;
for i = 1:iter
    subplot(iter,1,i)
    hold on;
    grid on;
    plot(tnl,xnl(i,:),'.','Markersize',10);
    plot(tl,xl(i,:),'linewidth',2);
    ylabel(labels{i},'Fontsize',14);
    if i==1
        legend('Measured', 'Predicted')
    end
end
xlabel('Time [s]','Fontsize',14);
suptitle(title)
set(gcf, 'Position', [100, 100, 1100, 730])

end