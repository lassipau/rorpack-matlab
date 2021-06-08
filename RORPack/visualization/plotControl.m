function plotControl(tgrid,CLsim,ContrSys,N,PrintFigureTitles)

% Plot the control input, if the arguments ContrSys and N are given
plot(tgrid,[zeros(size(ContrSys.K,1),N),ContrSys.K]*CLsim.xesol,'Linewidth',2);
set(gca,'xgrid','on','ygrid','on','tickdir','out','box','off')
if PrintFigureTitles == true
    title('Control input $u(t)$','Interpreter','latex','Fontsize',16)
end
set(gcf,'color',1/255*[252 247 255])
end