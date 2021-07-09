function PlotHeat2DSurfCase2(CLsim,spgrid,tgrid,plotind,zlims)
%
% Plot the solution of the 2D Heat equation
% CLsim = simulated state of the closed-loop system
% spgrid = spatial grid for the animation
% tgrid = t grid for the animation
% plotind = tgrid(plotind) is the time t at which the solution is plotted
% zlims = limits for the z-axis (optional)


xx = spgrid.xx;
yy = spgrid.yy;

% size of the state space X_N is N*M, here the grid has is N x M
N = size(xx,1);
M = size(xx,2);

zz = CLsim.xesol(1:N*M,plotind);

if max(max(abs(imag(zz)))) > 1e-8
  warning('Solution may contain imaginary parts that are ignored in the animation.')
end

zz = real(zz);

if nargin <= 4
  zlims(1) = min(min(min(zz)));
  zlims(2) = max(max(max(zz)));
end
axlims = [xx(1,1) xx(1,end) yy(1,1) yy(end,1) zlims];

% state at time tgrid(plotind)
surf(xx,yy,reshape(zz,M,N).');
axis(axlims)
caxis(zlims)
ylabel('$y$','Interpreter','latex','Fontsize',20)
xlabel('$x$','Interpreter','latex','Fontsize',20)
title(['The state of the system at $t=' num2str(round(tgrid(plotind),1)) '$.'],'interpreter', 'latex', 'fontsize', 20);
