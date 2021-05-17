function PlotLHMNCSurf(state,spgrid,tgrid,zlims)
% function PlotHeat2DSurf(state,spgrid)
%
% Plot the solution of the Timoshenko Beam, only the fourth component in
% the solution
% state = state of the closed-loop system at times 'tgrid'
% zlims = limits for the z-axis (optional)


% size of the state space X_N is N^2, here the grid has is N x N
N = length(spgrid)-1;

zz = state((3*N+1):(4*N),:);

if max(max(abs(imag(zz)))) > 1e-8
  warning('Solution may contain imaginary parts that are ignored in the animation.')
end

zz = real(zz);

if nargin <= 3
  zlims(1) = min(min(min(zz)));
  zlims(2) = max(max(max(zz)));
end
axlims = [tgrid(1) tgrid(end) spgrid(1) spgrid(end) zlims];
    
surf(tgrid,spgrid(2:end),zz);
axis(axlims)
caxis(zlims)
set(gca, 'YDir', 'reverse')

xlabel('$t$','Interpreter','latex','Fontsize',20)
ylabel('$\xi$','Interpreter','latex','Fontsize',20)

