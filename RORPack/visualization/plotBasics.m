function plotBasics(tgrid,yref,CLsim,ContrSys,N)

% Plot the output with the reference signal
if nargin > 3
    subplot(3,1,1)
else
    subplot(2,1,1)
end
hold off
cla
hold on
plot(tgrid,yref(tgrid),'Color',1.1*[0 0.447 0.741],'Linewidth',2);
plot(tgrid,CLsim.output,'Color', [0.85 0.325 0.098],'Linewidth',2);
title('Output $y(t)$ (red) and the reference $y_{ref}(t)$ (blue)','Interpreter','latex','Fontsize',16)
set(gca,'xgrid','off','tickdir','out','box','off')
if nargin > 3
    subplot(3,1,2)
else
    subplot(2,1,2)
end

% Plot the regulation error
plot(tgrid,CLsim.error,'Linewidth',2);
set(gca,'xgrid','on','ygrid','on','tickdir','out','box','off')
title('Regulation error $y(t)-y_{ref}(t)$','Interpreter','latex','Fontsize',16)
set(gcf,'color',1/255*[252 247 255])

% Plot the control input, if the arguments ContrSys and N are given
if nargin > 3
    subplot(3,1,3)
    plot(tgrid,[zeros(size(ContrSys.K,1),N),ContrSys.K]*CLsim.xesol,'Linewidth',2);
    set(gca,'xgrid','on','ygrid','on','tickdir','out','box','off')
    title('Control input $u(t)$','Interpreter','latex','Fontsize',16)
    set(gcf,'color',1/255*[252 247 255])
end