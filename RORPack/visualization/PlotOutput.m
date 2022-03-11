function PlotOutput(tgrid,yref,CLsim,PrintFigureTitles)

% Plot the output with the reference signal
hold off
cla
hold on
plot(tgrid,yref(tgrid),'Color',1.1*[0 0.447 0.741],'Linewidth',2);
plot(tgrid,CLsim.output,'Color', [0.85 0.325 0.098],'Linewidth',2);
if PrintFigureTitles == true
    title('Output $y(t)$ (red) and the reference $y_{ref}(t)$ (blue)','Interpreter','latex','Fontsize',16)
end
set(gca,'xgrid','off','tickdir','out','box','off')
end