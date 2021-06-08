function [MovAnim,zlims] = Anim1DHeat(state,spgrid,tgrid,BCtype,Tpause,record,zlims)
% function [MovAnim,zlims] = Anim1DHeat(state,spgrid,tgrid,BCtype,Tpause,record,zlims)
%
% Plot the solution (displacement w) of the controlled 1D Wave equation
% state = state of the closed-loop system at times 'tgrid'
% spgrid = spatial grid
% tgrid = grid for time
% BCtype = types of the boundary conditions, 'NN', 'ND', etc
% Tpause = pause time in the animation
% record = 1 if animation is recorded into a movie
% zlims = limits for the z-axis (optional)


if isequal(BCtype(1),'D')
  state = [zeros(1,size(state,2));state];
elseif ~isequal(BCtype(1),'N')
  error('Unknown boundary condition type')
end

if isequal(BCtype(2),'D')
  state = [state;zeros(1,size(state,2))];
elseif ~isequal(BCtype(2),'N')
  error('Unknown boundary condition type')
end

if max(max(abs(imag(state)))) > 1e-8
  warning('Solution may contain imaginary parts that are ignored in the plot.')
end

state = real(state);

if nargin <= 6
  zlims(1) = min(min(min(state)));
  zlims(2) = max(max(max(state)));
end
axlims = [spgrid(1) spgrid(end) zlims];
    


if record == 1
  
  MovAnim = struct('cdata',[],'colormap',[]);
  
  for ind = 1:size(state,2)
    
    plot(spgrid,state(:,ind).','Linewidth',2)
    axis(axlims)
    xlabel('$\xi$','Interpreter','latex','Fontsize',20)
    title(['Time $=\; ' num2str(tgrid(ind),'%.1f') '$'],'Interpreter','latex','Fontsize',20)
    set(gcf,'color',1/255*[252 247 255])
    drawnow
    MovAnim(ind) = getframe(gcf);
    pause(Tpause)
    
  end
  
else
  
  MovAnim = [];
  for ind = 1:size(state,2)
    
    plot(spgrid,state(:,ind).','Linewidth',2)
    axis(axlims)
    xlabel('$\xi$','Interpreter','latex','Fontsize',20)
    title(['Time $=\; ' num2str(tgrid(ind),'%.1f') '$'],'Interpreter','latex','Fontsize',20)
    set(gcf,'color',1/255*[252 247 255])
    drawnow
    pause(Tpause)
    
  end
end