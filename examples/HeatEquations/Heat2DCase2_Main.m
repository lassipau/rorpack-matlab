%% Robust control of a 2D heat equation with either the Observer Based
% Robust Controller (ObsBasedRC) or the Dual Observer Based Robust 
% Controller (DualObsBasedRC). 
% Neumann boundary control and Dirichlet boundary observation.
% Approximation with a Finite differences scheme.

% The system is unstable with a single unstable eigenvalue s=0.

% addpath(genpath('../RORPack/'))


cval = 1;
N = 15; 
M = 16;


% Initial state of the plant
%x0fun = @(x,y) zeros(size(x));
x0fun = @(x,y) 0.5*(1+cos(pi*(1-x))).*(1-1/4*cos(2*pi*y));
%x0fun = @(x,y) 0.5*(1+cos(pi*(1-x)));
%x0fun = @(x,y) cos(pi*(1-x));
%x0fun = @(x,y) 1/2*x.^2.*(3-2*x)-1;
%x0fun = @(x,y) 1/2*x.^2.*(3-2*x)-1/2;
%x0fun = @(x,y) 1*(1-x).^2.*(3-2*(1-x))-1;
%x0fun = @(x,y) .5*(1-x).^2.*(3-2*(1-x))-.5;
%x0fun = @(x,y) 1/4*(x.^3-1.5*x.^2)-1/4;
%x0fun = @(x,y) .2*x.^2.*(3-2*x)-.5;

[x0,spgrid,Sys] = ConstrHeat2DCase2(N,M,x0fun,cval);

% Define the reference and disturbance signals, and list the
% required frequencies in 'freqsReal'

%yref = @(t) sin(2*t)+.1*cos(6*t);
%yref = @(t) sin(2*t)+.2*cos(3*t);
%yref = @(t) ones(size(t));

%wdist = @(t) ones(size(t));
%wdist = @(t) sin(t);
%wdist = @(t) zeros(size(t));

% Case 1:
yref = @(t) sin(2*t)+0.1*cos(3*t);
%wdist = @(t) zeros(size(t));
%wdist = @(t) sin(2*t);

% Case 2:
% yref = @(t) ones(size(t));
% wdist = @(t) ones(size(t));

% Case 3:
% yref = @(t) sin(2*t)+.1*cos(6*t);
wdist = @(t) sin(t);

freqsReal = [0, 1, 2, 3];

% Check the consistency of the system definition
Sys = SysConsistent(Sys,yref,wdist,freqsReal);

%% Construct the controller 

% Observer-based robust controller
%
% Requires the construction of the stabilizing operators K21 and L to
% stabilize the single unstable eigenvalue s=0.
%
% 1 for interior points, 1/4 for corner points and 1/2 for other boundary
% points
K21 = -ones(1,N*M);
K21(1, 1:M) = K21(1, 1:M) / 2;
K21(1, (N-1)*M+1:end) = K21(1, (N-1)*M+1:end) / 2;
K21(1, 1:N:end) = K21(1, 1:N:end) / 2;
K21(1, M:N:end) = K21(1, M:N:end) / 2;
L = K21' * 10;
IMstabmarg = 0.5;
IMstabtype = 'LQR';
% IMstabtype = 'poleplacement';

ContrSys = ObserverBasedRC(freqsReal,Sys,K21,L,IMstabtype,IMstabmarg);


% % A Low-Gain 'Minimal' Robust Controller
% %
% % Even though the system is not impedance passive, the single unstable
% % eigenvalue s=0 can be stabilized with negative output feedback
% % u(t)=-kappa_S*y(t). This can be used in the design of the Low Gain Robust
% % Controller with a feedthrough operator -kappa_S*eye(dimY) (here dimY=1).
% 
% % Stabilizing output feedback gain
% kappa_S = -4;
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
%         Pval_nominal = Pappr(1i*freqsReal(ind));
%         Pvals{ind} = (eye(dimY)-Pval_nominal*Dc)\Pval_nominal;
%     end
% end
% 
% epsgain = [0.01,2];
% 
% % Construct the Low-Gain Robust Controller with a feedthrough term Dc.
% [ContrSys,epsgain] = LowGainRC(freqsReal,Pvals,epsgain,Sys,Dc);
% epsgain

%% Closed-loop construction and simulation

% Construct the closed-loop system
CLSys = ConstrCLSys(Sys,ContrSys);

% Print an approximate stability margin of the closed-loop system
stabmarg = CLStabMargin(CLSys)

% Define the initial state of the closed-loop system
% (the controller has zero initial state by default).
xe0 = [x0;zeros(size(ContrSys.G1,1),1)];

% Set the simulation length and define the plotting grid
Tend = 8;
tgrid = linspace(0,Tend,300);

% Simulate the closed-loop system
CLsim = SimCLSys(CLSys,xe0,yref,wdist,tgrid,[]);

%% Visualization

% Plot the (approximate) eigenvalues of the closed-loop system
figure(1)
% PlotEigs(CLSys.Ae,[-1 .3 -4 4]);
PlotEigs(CLSys.Ae,[-2 .3 -6 6]);

% Choose whether or not to print titles of the figures
PrintFigureTitles = true;

% Plot the controlled outputs, the tracking error norm, and 
% the control inputs
figure(2)
subplot(3,1,1)
PlotOutput(tgrid,yref,CLsim,PrintFigureTitles)
subplot(3,1,2)
PlotErrorNorm(tgrid,CLsim,PrintFigureTitles)
subplot(3,1,3)
PlotControl(tgrid,CLsim,PrintFigureTitles)
set(gcf,'color','white')

%% State of the controlled PDE

figure(3)
% PlotHeat2DSurf(x0,spgrid,[-1.4,1.4])
plotind = 155;
PlotHeat2DSurfCase2(CLsim,spgrid,tgrid,plotind)

%% Animation of the state of the controlled PDE

figure(4)
colormap jet
% No movie recording
[~,zlims] = AnimHeat2DCase2(CLsim,spgrid,tgrid,0.03,0);

% Movie recording
% [MovAnim,zlims] = AnimHeat2DCase2(CLsim,spgrid,tgrid,0,1);

%movie(MovAnim)

%% The reference signal

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
