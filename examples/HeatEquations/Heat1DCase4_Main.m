%% Heat equation on the interval [0,1] with 
% Dirichlet boundary control at x=1. Dirichlet boundary observation and
% Neumann boundary disturbance at x=0.
% Approximation with a Finite differences scheme.

% Case 4: Dirichlet boundary control at x=1, regulated output y(t) and a 
% Neumann boundary disturbance at x=0. The system is exponentially stable,
% but does not define a "wellposed" or "regular linear system" on the
% natural state space X=L^2(0,1), but instead a "Boundary Control System".
% The Low-Gain Robust Controller can be used due to the theory in the
% references Humaloja-Paunonen IEEE TAC 2018 and Humaloja-Kurula-Paunonen
% IEEE TAC 2019. Also other controllers can be constructed, but since the 
% current theory does not guarantee that the controller designs would work, 
% these are simulations are only for experimentation purposes. 
% That is, PROCEED WITH CAUTION! ;)

% addpath(genpath('../RORPack/'))

N = 50; 

% Initial state of the plant
%x0fun = @(x) zeros(size(x));
x0fun = @(x) 0.5*(1+cos(pi*(1-x)));
% x0fun = @(x) 0.5*(1+cos(pi*(1-x)));
%x0fun = @(x) 1/2*x.^2.*(3-2*x)-1;
%x0fun = @(x) 1/2*x.^2.*(3-2*x)-1/2;
%x0fun = @(x) 1*(1-x).^2.*(3-2*(1-x))-1;
%x0fun = @(x) .5*(1-x).^2.*(3-2*(1-x))-.5;
%x0fun = @(x) 1/4*(x.^3-1.5*x.^2)-1/4;
%x0fun = @(x) .2*x.^2.*(3-2*x)-.5;

% The spatially varying thermal diffusivity of the material
% cfun = @(t) ones(size(t));
% cfun = @(t) 1+t;
cfun = @(t) 1-2*t.*(1-2*t);
% cfun = @(t) 1+0.5*cos(5/2*pi*t);
% cfun = @(t) 0.3-0.6*t.*(1-t);

[x0,Sys,spgrid,BCtype] = ConstrHeat1DCase4(cfun,x0fun,N);

% Model = ss(Sys.A,Sys.B,Sys.C,Sys.D);
% tt=linspace(0,4*pi);
% utestfun = @(t) ones(size(t));
% utestfun = @(t) sin(t);
% [output,t,xx]=lsim(Model,utestfun(tt),tt);
% xx = [xx, utestfun(tt(:))]
% plot(spgrid,xx(end,:))
% surf(spgrid,tt,xx)
% 
% % plot(tt,(xx(:,2)-xx(:,1))*(N-1))

%yref = @(t) sin(2*t)+.1*cos(6*t);
%yref = @(t) sin(2*t)+.2*cos(3*t);
%yref = @(t) ones(size(t));

%wdist = @(t) ones(size(t));
%wdist = @(t) sin(t);
%wdist = @(t) zeros(size(t));

% Case 1:
yref = @(t) sin(2*t);%+.2*cos(3*t);
wdist = @(t) zeros(size(t));
% wdist = @(t) sin(2*t);

% Case 2:
% yref = @(t) ones(size(t));
% wdist = @(t) ones(size(t));

% Case 3:
% yref = @(t) sin(2*t)+.1*cos(6*t);
% wdist = @(t) sin(t);

freqsReal = [0, 2];

% Check the consistency of the system definition
Sys = SysConsistent(Sys,yref,wdist,freqsReal);

%% Construct the controller

% A Low-Gain 'Minimal' Robust Controller

Pappr = @(s) Sys.C*((s*eye(size(Sys.A,1))-Sys.A)\Sys.B)+Sys.D;
Pvals = cell(1,length(freqsReal));
for ind = 1:length(freqsReal)
  Pvals{ind} = Pappr(freqsReal(ind));
end

epsgain = [0.05,3];
% epsgain = 1;
[ContrSys,epsgain] = LowGainRC(freqsReal,Pvals,epsgain,Sys);
epsgain

% % An observer-based robust controller
% % DISCLAIMER: THE EXISTING THEORY DOES NOT GUARANTEE THAT THESE 
% % CONTROLLER WOULD WORK FOR THIS SYSTEM 
% % Stabilizing state feedback and output injection operators K and L
% % The plant is already stable, but stability margin can be improved.
% K = -7*[1, zeros(1,N-1)];
% %K = zeros(1,N);
% %K = -ones(1,N);
% %PlotEigs(full(Sys.A+Sys.B*K),[-20 1 -.3 .3])
% 
% L = -20*[zeros(N-1,1);2*(N-1)];
% % L = zeros(N,1);
% % PlotEigs(full(Sys.A+L*Sys.C),[-20 1 -.3 .3])
% 
% % ContrSys = ObserverBasedRC(freqs,Sys,K,L,'LQR',4);
% ContrSys = ObserverBasedRC(freqsReal,Sys,K,L,'poleplacement',4);
% % ContrSys = DualObserverBasedRC(freqs,Sys,K,L,'LQR',4);
% % ContrSys = DualObserverBasedRC(freqs,Sys,K,L,'poleplacement',4);

%% Closed-loop simulation and visualization of the results

% Construct the closed-loop system
CLSys = ConstrCLSys(Sys,ContrSys);

stabmarg = CLStabMargin(CLSys)

figure(1)
PlotEigs(CLSys.Ae,[-20 .3 -6 6]);


% Simulate the closed-loop system
xe0 = [x0;zeros(size(ContrSys.G1,1),1)];

Tend = 14;
tgrid = linspace(0,Tend,300);



CLsim = SimCLSys(CLSys,xe0,yref,wdist,tgrid,[]);

% Choose whther or not to print titles of the figures
PrintFigureTitles = true;

figure(2)
subplot(3,1,1)
PlotOutput(tgrid,yref,CLsim,PrintFigureTitles)
subplot(3,1,2)
PlotErrorNorm(tgrid,CLsim,PrintFigureTitles)
subplot(3,1,3)
PlotControl(tgrid,CLsim,PrintFigureTitles)

%%


% In plotting and animating the state,
% fill in the Dirichlet boundary condition x(1,t)=u(t) at x=1
spgrid_plot = [spgrid, 1];
state_plot = [CLsim.xesol(1:N,:);CLsim.control];
inputs = [zeros(size(ContrSys.K,1),N),ContrSys.K]*CLsim.xesol;
BCtype = 'NN';

figure(3)
colormap jet
Plot1DHeatSurf(state_plot,spgrid_plot,tgrid,BCtype)

%%
figure(4)
% No movie recording
[~,zlims] = Anim1DHeat(state_plot,spgrid_plot,tgrid,BCtype,0.03,0);

% Movie recording
% [MovAnim,zlims] = Anim1DHeat(state_plot,spgrid_plot,tgrid,BCtype,0.03,1);
%movie(MovAnim)

%%


figure(5)
tt = linspace(0,16,500);
plot(tt,yref(tt),'Color',1.1*[0 0.447 0.741],'Linewidth',3);
title('Reference signal $y_{ref}$','Interpreter','latex','Fontsize',16)
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
