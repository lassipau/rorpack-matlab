function Plot1DHeatSurf(state,spgrid,tgrid,BCtype,zlims)
% function Plot1DHeatSurf(state,spgrid,tgrid,zlims,BCtype)
%
% Plot the solution (temperature) of the controlled 1D heat equation
% state = state of the heat system at times 'tgrid'. One column for each
% time instance
% spgrid = spatial grid
% tgrid = grid for time
% BCtype = 'NN' for Neumann-Neumann, 'ND' for Neumann-Dirichlet etc
% zlims = limits for the z-axis (optional)

if isequal(BCtype(1),'D')
  state = [zeros(1,length(tgrid));state];
elseif ~isequal(BCtype(1),'N')
  error('Unknown boundary condition type')
end

if isequal(BCtype(2),'D')
  state = [state;zeros(1,length(tgrid))];
elseif ~isequal(BCtype(2),'N')
  error('Unknown boundary condition type')
end

if max(max(abs(imag(state)))) > 1e-8
  warning('Solution may contain imaginary parts that are ignored in the plot.')
end

state = real(state);

if nargin <= 5
  zlims(1) = min(min(min(state)));
  zlims(2) = max(max(max(state)));
end
axlims = [tgrid(1) tgrid(end) spgrid(1) spgrid(end) zlims];
    
surf(tgrid,spgrid,state);
axis(axlims)
caxis(zlims)
set(gca, 'YDir', 'reverse')

xlabel('$t$','Interpreter','latex','Fontsize',20)
ylabel('$\xi$','Interpreter','latex','Fontsize',20)

