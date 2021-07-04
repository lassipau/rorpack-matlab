%% Heat equation on the interval [0,1] with 
% Neumann boundary control and Dirichlet boundary observation 
% Approximation with a Finite differences scheme 

% Case 5: Similar to Case 1, but with two boundary inputs and outputs: 
% Neumann boundary control u_1(t) at x=0, and u_2(t) at x=1. Pointwise
% temperature measurements y_1(t) at x=0, and y_2(t) at x=1. Two input
% disturbances w_{dist,1}(t) at x=0 and w_{dist,2}(t) at x=1. The system is
% unstable (eigenvalue at 0), but is impedance passive and can be
% stabilized with negative output feedback.

addpath(genpath('../RORPack/'))

N = 70; 

% Initial state of the plant
%x0fun = @(x) zeros(size(x));
x0fun = @(x) 1*(1+cos(pi*(1-x)));
% x0fun = @(x) 3*(1-x)+x;
% x0fun = @(x) 0.5*(1+cos(pi*(1-x)));
%x0fun = @(x) 1/2*x.^2.*(3-2*x)-1;
%x0fun = @(x) 1/2*x.^2.*(3-2*x)-1/2;
%x0fun = @(x) 1*(1-x).^2.*(3-2*(1-x))-1;
%x0fun = @(x) .5*(1-x).^2.*(3-2*(1-x))-.5;
%x0fun = @(x) 1/4*(x.^3-1.5*x.^2)-1/4;
%x0fun = @(x) .2*x.^2.*(3-2*x)-.5;

% The spatially varying thermal diffusivity of the material
cfun = @(t) ones(size(t));
% cfun = @(t) 1+t;
% cfun = @(t) 1-2*t.*(1-2*t);
% cfun = @(t) 1+0.5*cos(5/2*pi*t);
% cfun = @(t) 0.3-0.6*t.*(1-t);

% Input profile for the distributed input disturbance (w_{dist,3}(t))
Bd_profile = @(xi) sin(pi/2*xi);
% Bd_profile = @(xi) zeros(size(xi));

[x0,Sys,spgrid,BCtype] = ConstrHeat1DCase5(cfun,x0fun,N,Bd_profile);

% Define the reference and disturbance signals: the system has two outputs
% and three disturbance inputs (two input disturbances and a distributed
% disturbance)
% Case 1:
yref = @(t) [sin(2*t);2*cos(3*t)];
wdist = @(t) [sin(6*t);-ones(size(t));0.2*sin(3*t)];
% wdist = @(t) zeros(3,size(t));

% Case 2:
% yref = @(t) [sin(2*t)+.3*cos(6*t);ones(size(t))];
% wdist = @(t) [sin(1*t);ones(size(t))+.3*sin(3*t+0.3);ones(size(t))];


freqsReal = [0, 1, 2, 3, 6];

% Check the consistency of the system definition
Sys = SysConsistent(Sys,yref,wdist,freqsReal);


%% Construct the controller

% % A Low-Gain 'Minimal' Robust Controller 
% %
% % Since the system is not impedance passive, the unstable eigenvalue s=0 
% % can be stabilized with negative output feedback u(t)=-kappa_S*y(t). This 
% % can be used in the design of the Low Gain Robust Controller with a 
% % feedthrough operator -kappa_S*eye(dimY) (here dimY=2).
% 
% % Stabilizing output feedback gain
% kappa_S = -1.5;
% % Controller feedthrough term 
% Dc = kappa_S*eye(size(Sys.C,1));
% 
% % Since the controller has a feedthrough term, the Pvals need to be
% % computed for the "prestabilized" system (A+B*Dc*C,B,C) (since D=0).
% % For nonzero frequencies we can use the formula 
% % P_{Dc}(i*w_k) = (I-P(i*w_k)*Dc)\P(i*w_k), but for s=0 we need to compute
% % this separately due to the fact that P(0) is not defined.
% 
% dimX = size(Sys.A,1);
% dimY = size(Sys.C,1);
% Pappr = @(s) Sys.C*((s*eye(dimX)-Sys.A)\Sys.B)+Sys.D;
% 
% Pvals = cell(1,length(freqsReal));
% for ind = 1:length(freqsReal)
%     if freqsReal(ind)==0
%         Pvals{ind} = Sys.C*((-Sys.A-Sys.B*Dc*Sys.C)\Sys.B)+Sys.D;
%     else
%         Pval_nominal = Pappr(freqsReal(ind));
%         Pvals{ind} = (eye(dimY)-Pval_nominal*Dc)\Pval_nominal;
%     end
% end
% 
% epsgain = [0.01,4];
% % epsgain = .5;
% [ContrSys,epsgain] = LowGainRC(freqsReal,Pvals,epsgain,Sys,Dc);
% epsgain



% A Passive Robust Controller
%
% Since the system is unstable, requires prestabilization with a
% controller feedthrough term 'Dc'.

% Negative feedback gain for output stabilization
kappa_S = -2.5;

dimY = size(Sys.C,1);
epsgain = [0.01,3];
% epsgain = .1;
[ContrSys,epsgain] = PassiveRC(freqsReal,dimY,epsgain,Sys,kappa_S*eye(dimY));
epsgain


% % An observer-based robust controller
% % Stabilizing state feedback and output injection operators K and L
% % These are chosen either based on collocated design or using LQR/LQG. 
% % Only the single unstable eigenvalue at s=0 needs to be stabilized.
% K = -2*Sys.C;
% % K = -lqr(Sys.A,Sys.B,0.1*eye(N),10*eye(2));
% PlotEigs(full(Sys.A+Sys.B*K),[-10 .1 -1 1])
% 
% L = -3*Sys.B;
% % L = -lqr(Sys.A',Sys.C',10*eye(N),eye(2))';
% % PlotEigs(full(Sys.A+L*Sys.C),[-20 1 -.3 .3])
% 
% ContrSys = ObserverBasedRC(freqsReal,Sys,K,L,'LQR',1);
% % ContrSys = ObserverBasedRC(freqsReal,Sys,K,L,'poleplacement',1);
% % ContrSys = DualObserverBasedRC(freqsReal,Sys,K,L,'LQR',1);
% % ContrSys = DualObserverBasedRC(freqsReal,Sys,K,L,'poleplacement',1);
 

%% Closed-loop simulation and visualization of the results

% Construct the closed-loop system
CLSys = ConstrCLSys(Sys,ContrSys);

stabmarg = CLStabMargin(CLSys)

figure(1)
PlotEigs(CLSys.Ae,[-5 .3 NaN NaN])


% Simulate the closed-loop system
xe0 = [x0;zeros(size(ContrSys.G1,1),1)];

Tend = 16;
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


figure(3)
colormap jet
Plot1DHeatSurf(CLsim.xesol(1:N,:),spgrid,tgrid,BCtype)

%%
figure(4)
% No movie recording
[~,zlims] = Anim1DHeat(CLsim.xesol(1:N,:),spgrid,tgrid,BCtype,0.03,0);

% Movie recording
% [MovAnim,zlims] = Anim1DHeat(CLsim.xesol(1:N,:),spgrid,tgrid,BCtype,0.03,1);

% movie(MovAnim)

%%


figure(5)
tt = linspace(0,16,500);
plot(tt,yref(tt),'Color',1.1*[0 0.447 0.741],'Linewidth',3);
title('Reference signal $y_{ref}$','Interpreter','latex','Fontsize',16)
set(gca,'xgrid','on','ygrid','on','tickdir','out','box','off')


%% Export movie to AVI

% AnimExport = VideoWriter('Heat1DCase3-animation.mp4','MPEG-4');
% AnimExport.Quality = 100;

% AnimExport = VideoWriter('Heat1DCase3-animation.avi','Uncompressed AVI');

% AnimExport.FrameRate = 15;
% open(AnimExport);
% writeVideo(AnimExport,MovAnim);
% close(AnimExport);
