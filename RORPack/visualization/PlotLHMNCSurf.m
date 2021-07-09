function PlotLHMNCSurf(profile,spgrid,tgrid,zlims)
% function PlotHeat2DSurf(state,spgrid)
%
% Plot the controlled deflection profile of the Timoshenko Beam. 
% profile =  the (precomputed) deflection profile w(\xi,t) of the beam at 
%            points in 'spgrid' (in rows) and at times 'tgrid' (in columns)
% spgrid  =  the spatial grid on [0,1]
% tgrid   =  the time grid
% zlims   =  limits for the z-axis (optional)


if max(max(abs(imag(profile)))) > 1e-8
  warning('Solution may contain imaginary parts that are ignored in the animation.')
end

profile = real(profile);

if nargin <= 3
  zlims(1) = min(min(min(profile)));
  zlims(2) = max(max(max(profile)));
end
axlims = [tgrid(1) tgrid(end) spgrid(1) spgrid(end) zlims];
    
surf(tgrid,spgrid,profile);
axis(axlims)
caxis(zlims)
set(gca, 'YDir', 'reverse')

xlabel('$t$','Interpreter','latex','Fontsize',20)
ylabel('$\xi$','Interpreter','latex','Fontsize',20)

