function [MovAnim,zlims] =...
    AnimLHMNC(profile,spgrid,tgrid,Tpause,record,zlims)
% function [MovAnim,zlims] =...
%   AnimLHMNC(profile,spgrid,tgrid,Tpause,record,zlims)
%
% Animate the controlled deflection profile of the Timoshenko Beam. 
% profile =  the (precomputed) deflection profile w(\xi,t) of the beam at 
%            points in 'spgrid' (in rows) and at times 'tgrid' (in columns)
% spgrid  =  the spatial grid on [0,1]
% tgrid   =  the time grid
% Tpause  =  pause time in the animation
% record  =  1 if animation is recorded into a movie
% zlims   =  limits for the z-axis (optional)

if max(max(abs(imag(profile)))) > 1e-8
  warning('Solution may contain imaginary parts that are ignored in the animation.')
end

profile = real(profile);

if nargin <= 5
  zlims(1) = min(min(min(profile)));
  zlims(2) = max(max(max(profile)));
end
axlims = [spgrid(1) spgrid(end) zlims];


if record == 1
  
  MovAnim = struct('cdata',[],'colormap',[]);
  
  for ind = 1:size(profile,2)
    
    plot(spgrid,profile(:,ind).','Linewidth',2)
    axis(axlims)
    xlabel('$\xi$','Interpreter','latex','Fontsize',20)
    title(['Time $=\; ' num2str(tgrid(ind),'%.1f') '$'],'Interpreter','latex','Fontsize',20)
    % set(gcf,'color',1/255*[252 247 255])
    drawnow
    MovAnim(ind) = getframe(gcf);
    pause(Tpause)
    
  end
  
else
  
  MovAnim = [];
  for ind = 1:size(profile,2)
    
    plot(spgrid,profile(:,ind).','Linewidth',2)
    axis(axlims)
    xlabel('$\xi$','Interpreter','latex','Fontsize',20)
    title(['Time $=\; ' num2str(tgrid(ind),'%.1f') '$'],'Interpreter','latex','Fontsize',20)
    % set(gcf,'color',1/255*[252 247 255])
    drawnow
    pause(Tpause)
    
  end
end