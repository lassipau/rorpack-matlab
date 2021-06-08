function [MovAnim,zlims] = Anim1DWaveSpectral(state,phin,spgrid,tgrid,Tpause,record,zlims)
% function [MovAnim,zlims] = Anim1DWaveSpectral(state,phin,spgrid,tgrid,zlims,Tpause,record)
%
% Plot the solution (displacement w) of the controlled 1D Wave equation
% state = state of the closed-loop system at times 'tgrid'
% phin = normalized basis of eigenfunctions of the Laplacian
% spgrid = spatial grid
% tgrid = grid for time
% Tpause = pause time in the animation
% record = 1 if animation is recorded into a movie
% zlims = limits for the z-axis (optional)


N = size(state,1)/2;
phinvals = phin(spgrid,(0:(N-1)).'); 
ww = state(1:2:end,:).'*phinvals;

if max(max(abs(imag(ww)))) > 1e-8
  warning('Solution may contain imaginary parts that are ignored in the plot.')
end

ww = real(ww);

if nargin <= 6
  zlims(1) = min(min(min(ww)));
  zlims(2) = max(max(max(ww)));
end

axlims = [spgrid(1) spgrid(end) zlims];
    
% surf(tgrid,spgrid,ww.');
% axis(axlims)
% caxis(zlims)
% set(gca, 'YDir', 'reverse')

% xlabel('$t$','Interpreter','latex','Fontsize',20)
% ylabel('$\xi$','Interpreter','latex','Fontsize',20)


if record == 1
  
  MovAnim = struct('cdata',[],'colormap',[]);
  
  for ind = 1:length(tgrid)
    
%     surf(xx,yy,reshape(zz(:,ind),N,N));
    plot(spgrid,ww(:,ind),'Linewidth',2)
    axis(axlims)
%     caxis(zlims)
    title(['Time $=\; ' num2str(tgrid(ind),'%.1f') '$'],'Interpreter','latex','Fontsize',20)
    xlabel('$\xi$','Interpreter','latex','Fontsize',20)
    set(gcf,'color',1/255*[252 247 255])
    drawnow
    MovAnim(ind) = getframe(gcf);
    pause(Tpause)
    
  end
  
else
  
  MovAnim = [];
  for ind = 1:size(ww,1)
    
%     surf(xx,yy,reshape(zz(:,ind),N,N));
    plot(spgrid,ww(ind,:),'Linewidth',2)
    axis(axlims)
%     caxis(zlims)
    title(['Time $=\; ' num2str(tgrid(ind),'%.1f') '$'],'Interpreter','latex','Fontsize',20)
    xlabel('$\xi$','Interpreter','latex','Fontsize',20)
    set(gcf,'color',1/255*[252 247 255])
    drawnow
    pause(Tpause)
    
  end
end