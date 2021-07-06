%% Robust control of a 2D heat equation on a rectangle with boundary control 
% The example is the simulation example from the article "Controller Design 
% for Robust Output Regulation of Regular Linear Systems" by L. Paunonen,
% IEEE TAC 2016.
% The system has two boundary inputs and two collocated boundary outputs,
% as well as an additional boundary disturbance input. The system is
% initially unstable, but its only unstable eigenvalue 0 is prestabilized 
% with negative output feedback A-B*C. The system is also impedance passive 
% (due to the collocated inputs and outputs).

% addpath(genpath('../RORPack/'))

% Parameters for this example.
N = 31; 

% Initial state of the plant
x0fun = @(x,y) zeros(size(x));
%x0fun = @(x,y) 0.5*(1+cos(pi*(1-x))).*(1-1/4*cos(2*pi*y));
% x0fun = @(x,y) 0.5*(1+cos(pi*(1-x)));
%x0fun = @(x,y) 1/2*x.^2.*(3-2*x)-1;
%x0fun = @(x,y) 1/2*x.^2.*(3-2*x)-1/2;
%x0fun = @(x,y) 1*(1-x).^2.*(3-2*(1-x))-1;
%x0fun = @(x,y) .5*(1-x).^2.*(3-2*(1-x))-.5;
%x0fun = @(x,y) 1/4*(x.^3-1.5*x.^2)-1/4;
%x0fun = @(x,y) .2*x.^2.*(3-2*x)-.5;

% Construct the system.
[x0,spgrid,Sys] = ConstrHeat2DCase1(1,x0fun,N);

% Define the reference and disturbance signals, and list the
% required frequencies in 'freqsReal'

% Case 1:
yref = @(t) [(-1) * ones(size(t));cos(pi*t)];
wdist = @(t) zeros(size(t));
freqsReal = [0, pi];

% Case 2:
% yref = @(t) ones(size(t));
% wdist = @(t) ones(size(t));
% freqsReal = 0;

% Case 3:
% yref = @(t) sin(2*t)+.1*cos(6*t);
% wdist = @(t) sin(t);
% freqsReal = [0,1,2,6];

% Check the consistency of the system definition
Sys = SysConsistent(Sys,yref,wdist,freqsReal);

%% Construct the controller

% A Low-Gain 'Minimal' Robust Controller
%
Pappr = @(s) Sys.C*((s*eye(size(Sys.A,1))-Sys.A)\Sys.B)+Sys.D;

Pvals = cell(1,length(freqsReal));
for ind = 1:length(freqsReal)
  Pvals{ind} = Pappr(freqsReal(ind));
end

epsgain = [0.01,4];
% epsgain = .5;
[ContrSys,epsgain] = LowGainRC(freqsReal,Pvals,epsgain,Sys);
epsgain

%% Closed-loop construction and simulation

% Construct the closed-loop system
CLSys = ConstrCLSys(Sys,ContrSys);

% Print an approximate stability margin of the closed-loop system
stabmarg = CLStabMargin(CLSys)

% Define the initial state of the closed-loop system
% (the controller has zero initial state by default).
xe0 = [x0;zeros(size(ContrSys.G1,1),1)];

% Set the simulation length and define the plotting grid
Tend = 16;
tgrid = linspace(0,Tend,300);

% Simulate the closed-loop system
CLsim = SimCLSys(CLSys,xe0,yref,wdist,tgrid,[]);

%% Visualization

% Plot the (approximate) eigenvalues of the closed-loop system
figure(1)
PlotEigs(CLSys.Ae,[-1 .3 -4 4]);
% PlotEigs(CLSys.Ae,[-2 .3 -6 6]);

% Choose whther or not to print titles of the figures
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


%% State of the controlled PDE

% figure(3)
% colormap jet
% % PlotHeat2DSurf(x0,spgrid,[-1.4,1.4])
% PlotHeat2DSurf(x0,spgrid,zlims)

%% Animation of the state of the controlled PDE

figure(4)
colormap jet
% No movie recording
[~,zlims] = AnimHeat2DCase1(CLsim,spgrid,tgrid,0.03,0);

% Movie recording
% [MovAnim,zlims] = AnimHeat2DCase1(CLsim,spgrid,tgrid,0,1);

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
