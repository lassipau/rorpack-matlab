function PlotHeat2DSurfCase1(CLsim,spgrid,tgrid)
%
% Plot the solution of the 2D Heat equation
% 
% 
% 

xx = spgrid.xx;
yy = spgrid.yy;

% size of the state space X_N is N^2, here the grid has is N x N
N = size(xx,1);

xesol = NormalizeHeat2DData(N, deval(CLsim.solstruct,tgrid), xx, yy);

zz = xesol(1:N^2,:);

if max(max(abs(imag(zz)))) > 1e-8
  warning('Solution may contain imaginary parts that are ignored in the animation.')
end

zz = real(zz);

zlims(1) = min(min(min(zz)));
zlims(2) = max(max(max(zz)));
axlims = [xx(1,1) xx(1,end) yy(1,1) yy(end,1) zlims];

% state at time ind
ind = 150;
surf(xx,yy,reshape(zz(:,ind),N,N));
axis(axlims)
caxis(zlims)
ylabel('$y$','Interpreter','latex','Fontsize',20)
xlabel('$x$','Interpreter','latex','Fontsize',20)
title(['The state of the system at $t=' num2str(round(tgrid(ind),1)) '$.'],'interpreter', 'latex', 'fontsize', 20);
