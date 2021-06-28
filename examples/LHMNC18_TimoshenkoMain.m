%% Timoshenko beam for the LHMNC18 presentation
% Approximation with a Finite differences scheme 

% Stable and passive system, collocated I/O

addpath(genpath('../RORPack/'))




N = 50; 

% Initial state of the plant
x0fun = @(x) zeros(size(x));
%x0fun = @(x,y) 0.5*(1+cos(pi*(1-x))).*(1-1/4*cos(2*pi*y));
% x0fun = @(x,y) 0.5*(1+cos(pi*(1-x)));
%x0fun = @(x,y) 1/2*x.^2.*(3-2*x)-1;
%x0fun = @(x,y) 1/2*x.^2.*(3-2*x)-1/2;
%x0fun = @(x,y) 1*(1-x).^2.*(3-2*(1-x))-1;
%x0fun = @(x,y) .5*(1-x).^2.*(3-2*(1-x))-.5;
%x0fun = @(x,y) 1/4*(x.^3-1.5*x.^2)-1/4;
%x0fun = @(x,y) .2*x.^2.*(3-2*x)-.5;

[x0,spgrid,Sys] = ConstrLHMNC18(x0fun,N);

%yref = @(t) sin(2*t)+.1*cos(6*t);
%yref = @(t) sin(2*t)+.2*cos(3*t);
%yref = @(t) ones(size(t));

%wdist = @(t) ones(size(t));
%wdist = @(t) sin(t);
%wdist = @(t) zeros(size(t));

% Case 1:
yref = @(t) sin(2*t)+.5*cos(1*t);
wdist = @(t) zeros(size(t));

% Case 2:
% yref = @(t) ones(size(t));
% wdist = @(t) ones(size(t));

% Case 3:
% yref = @(t) sin(2*t)+.1*cos(6*t);
% wdist = @(t) sin(t);


% freqs = [-3i -2i -1i 0 1i 2i 3i];
freqsReal = [1 2];
freqs = 1i*freqsReal;

% dimX = size(Sys.A,1);
% Pappr = @(s) Sys.C*((s*eye(dimX)-Sys.A)\Sys.B)+Sys.D;
% 
% Pvals = cell(1,length(freqs));
% for ind = 1:length(freqs)
%   Pvals{ind} = Pappr(freqs(ind));
% end
% 
% epsgainrange = [10,50];
% epsgain = 13;
%[ContrSys,epsgain] = LowGainRC(freqs,Pvals,epsgain,Sys);
% [ContrSys,epsgain] = PassiveRC(freqsReal,Pvals,epsgainrange,Sys);
% [ContrSys,epsgain] = PassiveRC(freqsReal,Pvals,epsgain,Sys);

% epsgain

% Choose stabilizing operators K and L

% An observer-based robust controller
% Stabilizing state feedback and output injection operators K and L
% These are chosen based on collocated design. 
K = -.5*Sys.B';
% PlotEigs(full(Sys.A+Sys.B*K),[-1 .1 NaN NaN])

L = -.5*Sys.C';
% PlotEigs(full(Sys.A+L*Sys.C))

% ContrSys = ObserverBasedRC(freqsReal,Sys,K,L,'poleplacement',0.5);
ContrSys = ObserverBasedRC(freqsReal,Sys,K,L,'LQR',0.5);


CLSys = ConstrCLSys(Sys,ContrSys);
PlotEigs(full(CLSys.Ae))

stabmarg = CLStabMargin(CLSys)

x0 = zeros(4*N,1);
xe0 = [x0;zeros(size(ContrSys.G1,1),1)];

Tend = 16;
tgrid = linspace(0,Tend,300);



CLsim = SimCLSys(CLSys,xe0,yref,wdist,tgrid,[]);

figure(1)
subplot(2,1,1)
hold off
cla
hold on
plot(tgrid,yref(tgrid),'Color',1.1*[0 0.447 0.741],'Linewidth',2);
plot(tgrid,CLsim.output,'Color', [0.85 0.325 0.098],'Linewidth',2);
title('Output $y(t)$ (red) and the reference $y_{ref}(t)$ (blue)','Interpreter','latex','Fontsize',16)
set(gca,'xgrid','off','tickdir','out','box','off')
subplot(2,1,2)
plot(tgrid,CLsim.error,'Linewidth',2);
set(gca,'xgrid','on','ygrid','on','tickdir','out','box','off')
title('Tracking  error $y(t)-y_{ref}(t)$','Interpreter','latex','Fontsize',16)
%set(gcf,'color',1/255*[252 247 255])


figure(3)
colormap jet
%PlotHeat2DSurf(x0,spgrid,[-1.4,1.4])
PlotLHMNCSurf(CLsim.xesol,spgrid,tgrid,[-9 9])
%PlotLHMNCSurf(CLsim.xesol(205:end,:),spgrid,tgrid,[-9 9])
set(gca,'ztick',-8:4:8);
%colormap jet

% figure(4)
% tt = linspace(0,16,500)
% plot(tt,yref(tt),'Color',1.1*[0 0.447 0.741],'Linewidth',3);
% set(gca,'xgrid','on','ygrid','on','tickdir','out','box','off')

%%


% figure(2)
% colormap jet
% No movie recording
% [~,zlims] = AnimHeat2Dtest1(CLsim,spgrid,tgrid,0.03,0);

% Movie recording
% [MovAnim,zlims] = AnimHeat2Dtest1(CLsim,spgrid,tgrid,0,1);

%movie(MovAnim)

figure(3)
colormap jet
%PlotHeat2DSurf(x0,spgrid,[-1.4,1.4])
PlotHeat2DSurf(x0,spgrid,zlims)

figure(4)
tt = linspace(0,16,500);
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
% AnimExport = VideoWriter('Case3-animation.avi','Uncompressed AVI');

% AnimExport.FrameRate = 15;
% open(AnimExport);
% writeVideo(AnimExport,MovAnim);
% close(AnimExport);
