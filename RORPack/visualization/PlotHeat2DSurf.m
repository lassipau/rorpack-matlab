function PlotHeat2DSurf(state,spgrid,zlims)
% function PlotHeat2DSurf(state,spgrid)
%
% Plot the solution of the 2D Heat equation
% state = state of the closed-loop system 
% zlims = limits for the z-axis (optional)


xx = spgrid.xx;
yy = spgrid.yy;

% size of the state space X_N is N^2, here the grid has is N x N
N = size(xx,1);

zz = state(1:N^2);

if max(max(abs(imag(zz)))) > 1e-8
  warning('Solution may contain imaginary parts that are ignored in the animation.')
end

zz = real(zz);

if nargin <= 2
  zlims(1) = min(min(min(zz)));
  zlims(2) = max(max(max(zz)));
end
axlims = [xx(1,1) xx(1,end) yy(1,1) yy(end,1) zlims];
    
h = surf(xx,yy,reshape(zz,N,N));
% set(h,'linewidth',0.5,'edgecolor',0.5*[1 1 1])
axis(axlims)
caxis(zlims)

ylabel('$y$','Interpreter','latex','Fontsize',20)

