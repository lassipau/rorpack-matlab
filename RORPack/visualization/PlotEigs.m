function PlotEigs(A,axlim)
% function PlotEigs(A,axlims)
%
% Plots the eigenvalues of A
% If 'axlim' is not given, limits determined from the spectrum. If some
% components of 'axlim' are set as NaN, then also those limits will also be 
% determined from the spectrum.

Aspec = eig(full(A));

if nargin == 1
  relim = [min(real(Aspec)) max(real(Aspec))];
  imlim = [min(imag(Aspec)) max(imag(Aspec))];
  axlim = [relim imlim];
else
  if isnan(axlim(1))
    axlim(1) = min(real(Aspec));
  end
  if isnan(axlim(2))
    axlim(2) = max(real(Aspec));
  end
  if isnan(axlim(3))
    axlim(3) = min(imag(Aspec));
  end
  if isnan(axlim(4))
    axlim(4) = max(imag(Aspec));
  end
end
  

hold off
cla
hold on
plot(real(Aspec),imag(Aspec),'r.','Markersize',15)
% set the limits of the plot
axis(axlim)

% plot the axes
plot(axlim(1:2),[0 0],'k',[0 0],axlim(3:4),'k','Linewidth',1)


maxreal = num2str(max(real(Aspec)));
maxrealLF = num2str(max(real(Aspec(abs(imag(Aspec))<100))));
title(['Largest real part = $' maxreal '$ (Low freq = $ ' maxrealLF '$)' ],'Interpreter','Latex')