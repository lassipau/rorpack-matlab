%% Heat equation on a square with Neumann boundary conditions and 
% Neumann boundary control and Dirichlet boundary observation 
% Approximation with a Finite differences scheme 

% Unstable system, stabilization by stabilizing the only unstable
% eigenvalue =0

addpath(genpath('../RORPack/'))

N = 16; 

% Initial state of the plant
%x0fun = @(x,y) zeros(size(x));
x0fun = @(x,y) 0.5*(1+cos(pi*(1-x))).*(1-1/4*cos(2*pi*y));
% x0fun = @(x,y) 0.5*(1+cos(pi*(1-x)));
%x0fun = @(x,y) 1/2*x.^2.*(3-2*x)-1;
%x0fun = @(x,y) 1/2*x.^2.*(3-2*x)-1/2;
%x0fun = @(x,y) 1*(1-x).^2.*(3-2*(1-x))-1;
%x0fun = @(x,y) .5*(1-x).^2.*(3-2*(1-x))-.5;
%x0fun = @(x,y) 1/4*(x.^3-1.5*x.^2)-1/4;
%x0fun = @(x,y) .2*x.^2.*(3-2*x)-.5;

[x0,spgrid,Sys] = ConstrHeat2Dtest1(1,x0fun,N);

%yref = @(t) sin(2*t)+.1*cos(6*t);
%yref = @(t) sin(2*t)+.2*cos(3*t);
%yref = @(t) ones(size(t));

%wdist = @(t) ones(size(t));
%wdist = @(t) sin(t);
%wdist = @(t) zeros(size(t));

% Case 1:
yref = @(t) sin(2*t);%+.2*cos(3*t);
wdist = @(t) zeros(size(t));
wdist = @(t) sin(2*t);

% Case 2:
% yref = @(t) ones(size(t));
% wdist = @(t) ones(size(t));

% Case 3:
% yref = @(t) sin(2*t)+.1*cos(6*t);
% wdist = @(t) sin(t);


freqs = [-3i -2i -1i 0 1i 2i 3i];

if max(abs(real(freqs)))>0 && max(abs(imag(freqs)))>0
  error('nonzero real parts in frequencies!')
elseif max(abs(imag(freqs)))>0
  freqsReal = unique(abs(freqs));
end

freqsReal = [0 1 2 3 6];
freqs = 1i*freqsReal;

dimX = size(Sys.A,1);
Pappr = @(s) Sys.C*((s*eye(dimX)-Sys.A)\Sys.B)+Sys.D;

Pvals = cell(1,length(freqs));
for ind = 1:length(freqs)
  Pvals{ind} = Pappr(freqs(ind));
end

epsgainrange = [0.01,4];
% epsgain = .1;
%[ContrSys,epsgain] = ConstrContrLG(freqs,Pvals,epsgain,Sys);

[ContrSys,epsgain] = ConstrContrLGReal(freqsReal,Pvals,epsgainrange,Sys);
epsgain

% ContrSys = ConstrContrObsBasedReal(freqsReal,Sys);

CLSys = ConstrCLSys(Sys,ContrSys);

stabmarg = CLStabMargin(CLSys)
figure(2)
PlotEigs(CLSys.Ae,[-1 .3 -4 4]);

% PlotEigs(CLSys.Ae,[-2 .3 -6 6]);
%%


xe0 = [x0;zeros(size(ContrSys.G1,1),1)];

Tend = 30;
tgrid = linspace(0,Tend,300);



CLsim = SimCLSys(CLSys,xe0,yref,wdist,tgrid,[]);

figure(1)
subplot(3,1,1)
hold off
cla
hold on
plot(tgrid,yref(tgrid),'Color',1.1*[0 0.447 0.741],'Linewidth',2);
plot(tgrid,CLsim.output,'Color', [0.85 0.325 0.098],'Linewidth',2);
title('Output $y(t)$ (red) and the reference $y_{ref}(t)$ (blue)','Interpreter','latex','Fontsize',16)
set(gca,'xgrid','off','tickdir','out','box','off')
subplot(3,1,2)
plot(tgrid,CLsim.error,'Linewidth',2);
set(gca,'xgrid','on','ygrid','on','tickdir','out','box','off')
title('Regulation error $y(t)-y_{ref}(t)$','Interpreter','latex','Fontsize',16)
set(gcf,'color',1/255*[252 247 255])
subplot(3,1,2)
plot(tgrid,CLsim.error,'Linewidth',2);
set(gca,'xgrid','on','ygrid','on','tickdir','out','box','off')
title('Regulation error $y(t)-y_{ref}(t)$','Interpreter','latex','Fontsize',16)
set(gcf,'color',1/255*[252 247 255])
subplot(3,1,3)
% Plot the control input
plot(tgrid,[zeros(size(ContrSys.K,1),N^2),ContrSys.K]*CLsim.xesol,'Linewidth',2);
set(gca,'xgrid','on','ygrid','on','tickdir','out','box','off')
title('Control input $u(t)$','Interpreter','latex','Fontsize',16)
set(gcf,'color',1/255*[252 247 255])


%%


figure(2)
colormap jet
% No movie recording
% [~,zlims] = AnimHeat2Dtest1(CLsim,spgrid,tgrid,0.03,0);

% Movie recording
[MovAnim,zlims] = AnimHeat2Dtest1(CLsim,spgrid,tgrid,0,1);

%movie(MovAnim)

%%

figure(3)
colormap jet
%PlotHeat2DSurf(x0,spgrid,[-1.4,1.4])
PlotHeat2DSurf(x0,spgrid,zlims)

figure(4)
tt = linspace(0,16,500)
plot(tt,yref(tt),'Color',1.1*[0 0.447 0.741],'Linewidth',3);
set(gca,'xgrid','on','ygrid','on','tickdir','out','box','off')


%% Export movie to AVI

%AnimExport = VideoWriter('heat2Danim.avi','Uncompressed AVI');
% AnimExport = VideoWriter('Case1-animation.mp4','MPEG-4');
% AnimExport = VideoWriter('Case2-animation.mp4','MPEG-4');
% AnimExport = VideoWriter('Case3-animation.mp4','MPEG-4');
% AnimExport.Quality = 100;

% AnimExport = VideoWriter('Case1-animation.avi','Uncompressed AVI');
% AnimExport = VideoWriter('Case2-animation.avi','Uncompressed AVI');
AnimExport = VideoWriter('Case3-animation.avi','Uncompressed AVI');

AnimExport.FrameRate = 15;
open(AnimExport);
writeVideo(AnimExport,MovAnim);
close(AnimExport);
