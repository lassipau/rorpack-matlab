function Plot1DWaveSurf(state,phin,spgrid,tgrid,zlims)
% function PlotHeat2DSurf(state,spgrid)
%
% Plot the solution (displacement w) of the controlled 1D Wave equation
% state = state of the closed-loop system at times 'tgrid'
% phin = normalized basis of eigenfunctions of the Laplacian
% spgrid = spatial grid
% tgrid = grid for time
% zlims = limits for the z-axis (optional)


N = size(state,1)/2;
phinvals = phin(spgrid,(0:(N-1)).'); 
ww = state(1:2:end,:).'*phinvals;

if max(max(abs(imag(ww)))) > 1e-8
  warning('Solution may contain imaginary parts that are ignored in the plot.')
end

ww = real(ww);

if nargin <= 4
  zlims(1) = min(min(min(ww)));
  zlims(2) = max(max(max(ww)));
end
axlims = [tgrid(1) tgrid(end) spgrid(1) spgrid(end) zlims];
    
surf(tgrid,spgrid,ww.');
axis(axlims)
caxis(zlims)
set(gca, 'YDir', 'reverse')

xlabel('$t$','Interpreter','latex','Fontsize',20)
ylabel('$\xi$','Interpreter','latex','Fontsize',20)

