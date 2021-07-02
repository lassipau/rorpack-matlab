function [MovAnim,zlims] = AnimHeat2DCase1(CLsim,spgrid,tgrid,Tpause,record)
% function AnimHeat2Dtest1(CLsim,tgrid)
%
% Animate the solution of the 2D Heat equation
% CLsim = simulated state of the closed-loop system 
% tgrid = t grid for the animation
% Tpause = length of the pause between frames
% record = toggle movie recording on (1) or off (0)
% zlims = limits for the z-axis


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


if record == 1
  
  MovAnim = struct('cdata',[],'colormap',[]);
  
  for ind = 1:length(tgrid)
    
    surf(xx,yy,reshape(zz(:,ind),N,N));
    axis(axlims)
    caxis(zlims)
    xlabel('$x$','Interpreter','latex','Fontsize',20)
    ylabel('$y$','Interpreter','latex','Fontsize',20)
    title(['Time $=\; ' num2str(tgrid(ind),'%.1f') '$'],'Interpreter','latex','Fontsize',20)
    set(gcf,'color',1/255*[252 247 255])
    drawnow
    MovAnim(ind) = getframe(gcf);
    pause(Tpause)
    
  end
  
else
  
  MovAnim = [];
  for ind = 1:length(tgrid)
    
    surf(xx,yy,reshape(zz(:,ind),N,N));
    axis(axlims)
    caxis(zlims)
    xlabel('$x$','Interpreter','latex','Fontsize',20)
    ylabel('$y$','Interpreter','latex','Fontsize',20)
    title(['Time $=\; ' num2str(tgrid(ind),'%.1f') '$'],'Interpreter','latex','Fontsize',20)
    set(gcf,'color',1/255*[252 247 255])
    drawnow
    pause(Tpause)
    
  end
end



