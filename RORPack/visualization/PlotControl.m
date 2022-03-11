function PlotControl(tgrid,CLsim,PrintFigureTitles)

% Plot the control input, this is saved in the CLsim structure
plot(tgrid,CLsim.control,'Linewidth',2);
set(gca,'xgrid','on','ygrid','on','tickdir','out','box','off')
if PrintFigureTitles == true
    title('Control input $u(t)$','Interpreter','latex','Fontsize',16)
end
set(gcf,'color',1/255*[252 247 255])
end