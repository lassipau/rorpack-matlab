function plotErrorNorm(tgrid,CLsim,PrintFigureTitles)

% Plot the regulation error
hold off
cla
hold on
plot(tgrid,CLsim.error,'Linewidth',2);
set(gca,'xgrid','on','ygrid','on','tickdir','out','box','off')
if PrintFigureTitles == true
    title('Regulation error $y(t)-y_{ref}(t)$','Interpreter','latex','Fontsize',16)
end
set(gcf,'color',1/255*[252 247 255])