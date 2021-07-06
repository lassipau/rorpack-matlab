%% Heat equation on the interval [0,1] with 
% Neumann boundary control and Dirichlet boundary observation 
% Approximation with a Finite differences scheme 

% Case 2: Neumann boundary control at x=0, regulated output y(t) and a 
% Dirichlet boundary condition at x=1
% Unstable system, stabilization by stabilizing the only unstable
% eigenvalue =0

addpath(genpath('../RORPack/'))

N = 50;

% Initial state of the plant
%x0fun = @(x) zeros(size(x));
x0fun = @(x) 0.5*(1+cos(pi*(1-x)));
%x0fun = @(x) 3*(1-x)+x;
%x0fun = @(x) 0.5*(1+cos(pi*(1-x)));
%x0fun = @(x) 1/2*x.^2.*(3-2*x)-1;
%x0fun = @(x) 1/2*x.^2.*(3-2*x)-1/2;
%x0fun = @(x) 1*(1-x).^2.*(3-2*(1-x))-1;
%x0fun = @(x) .5*(1-x).^2.*(3-2*(1-x))-.5;
%x0fun = @(x) 1/4*(x.^3-1.5*x.^2)-1/4;
%x0fun = @(x) .2*x.^2.*(3-2*x)-.5;

% The spatially varying thermal diffusivity of the material
% cfun = @(t) ones(size(t));
% cfun = @(t) 1+t;
% cfun = @(t) 1-2*t.*(1-2*t);
cfun = @(t) 1+0.5*cos(5/2*pi*t);
% cfun = @(t) 0.3-0.6*t.*(1-t);

[x0,Sys,spgrid,BCtype] = ConstrHeat1DCase2(cfun,x0fun,N);

% Model = ss(Sys.A,Sys.B,Sys.C,Sys.D);
% tt=linspace(0,4);
% [output,t,xx]=lsim(Model,ones(size(tt)),tt);
% % plot(tt,(xx(:,2)-xx(:,1))*(N-1))
% surf(linspace(0,1,50),tt,xx)

% Define the reference and disturbance signals, and list the
% required frequencies in 'freqsReal'

%yref = @(t) sin(2*t)+.1*cos(6*t);
%yref = @(t) sin(2*t)+.2*cos(3*t);
%yref = @(t) ones(size(t));

%wdist = @(t) ones(size(t));
%wdist = @(t) sin(t);
%wdist = @(t) zeros(size(t));

% Case 1:
yref = @(t) sin(2*t) + .5*cos(3*t);
wdist = @(t) zeros(size(t));
% wdist = @(t) sin(6*t);

% Case 2:
% yref = @(t) ones(size(t));
% wdist = @(t) ones(size(t));

% Case 3:
% yref = @(t) sin(2*t)+.1*cos(6*t);
% wdist = @(t) sin(t);


% freqs = [-3i -2i -1i 0 1i 2i 3i];
freqsReal = [0 1 2 3 6];

% Sys.A = Sys.A+2*pi^2*Sys.B*Sys.Cm;
% PlotEigs(full(Sys.A),[-20 1 -.3 .3])

% eig(full(Sys.A))

% Check the consistency of the system definition
Sys = SysConsistent(Sys,yref,wdist,freqsReal);

%% Construct the controller

% A Low-Gain 'Minimal' Robust Controller
% 
% Compute the required transfer function values P(iw_k) using Chebfun: This
% can be done by solving a boundary value problem related to each
% frequency, and evaluating the solution at x=0.
Pvals = cell(1,length(freqsReal));
for ind = 1:length(freqsReal)
    s = freqsReal(ind);
    cb_A = chebop(0,1);
    cb_A.op = @(x,w) 1i*s*w-diff(cfun(x)*diff(w));
    cb_A.lbc = @(w) diff(w)+1;
    cb_A.rbc = @(w) w;
    w = cb_A\0;
    
    Pvals{ind} = w(0);
end

% % Alternative: Approximation using the finite difference approximation (if
% % Chebfun is not available)
% Pappr = @(s) Sys.C*((s*eye(size(Sys.A,1))-Sys.A)\Sys.B)+Sys.D;
% Pvals = cell(1,length(freqsReal));
% for ind = 1:length(freqsReal)
%    Pvals{ind} = Pappr(1i*freqsReal(ind));
% end
epsgainrange = [0.01,3];
% epsgain = .1;
[ContrSys,epsgain] = LowGainRC(freqsReal,Pvals,epsgainrange,Sys);
epsgain

% A Passive Robust Controller

% Compute the required transfer function values P(iw_k) using Chebfun: This
% can be done by solving a boundary value problem related to each
% frequency, and evaluating the solution at x=0.
Pvals = cell(1,length(freqsReal));
for ind = 1:length(freqsReal)
    s = freqsReal(ind);
    cb_A = chebop(0,1);
    cb_A.op = @(x,w) 1i*s*w-diff(cfun(x)*diff(w));
    cb_A.lbc = @(w) diff(w)+1;
    cb_A.rbc = @(w) w;
    w = cb_A\0;
    
    Pvals{ind} = w(0);
end
Pappr = @(s) Sys.C*((s*eye(size(Sys.A,1))-Sys.A)\Sys.B)+Sys.D;
Pvals = cell(1,length(freqsReal));
for ind = 1:length(freqsReal)
   Pvals{ind} = Pappr(1i*freqsReal(ind));
end

epsgainrange = [0.01,3];
% epsgain = .1;
[ContrSys,epsgain] = PassiveRC(freqsReal,Pvals,epsgainrange,Sys);
epsgain

% An observer-based robust controller or
% a dual observer-based robust controller
% Stabilizing state feedback and output injection operators K and L
% These are chosen based on collocated design. Only the single unstable
% eigenvalue at s=0 needs to be stabilized
% K = zeros(1,N);
% K = 1/(N-1)*[0.5, ones(1,N-2), 0.5];
% PlotEigs(full(Sys.A+Sys.B*K),[NaN .1 -.3 .3])

% L = -ones(N,1);
% PlotEigs(full(Sys.A+L*Sys.C),[-20 1 -.3 .3])

% ContrSys = ObserverBasedRC(freqsReal,Sys,K,L,'LQR',0.5);
% ContrSys = ObserverBasedRC(freqsReal,Sys,K,L,'poleplacement',0.5);
% ContrSys = DualObserverBasedRC(freqsReal,Sys,K,L,'LQR',0.5);
% ContrSys = DualObserverBasedRC(freqsReal,Sys,K,L,'poleplacement',0.5);

%% Closed-loop construction and simulation


CLSys = ConstrCLSys(Sys,ContrSys);

stabmarg = CLStabMargin(CLSys)

figure(1)
PlotEigs(CLSys.Ae,[-20 .3 -6 6])

xe0 = [x0;zeros(size(ContrSys.G1,1),1)];

Tend = 14;
tgrid = linspace(0,Tend,300);


CLsim = SimCLSys(CLSys,xe0,yref,wdist,tgrid,[]);

%% Visualization

% Choose whther or not to print titles of the figures
PrintFigureTitles = true;

% Plot the output, reference error and control
figure(2)
subplot(3,1,1)
PlotOutput(tgrid,yref,CLsim,PrintFigureTitles)
subplot(3,1,2)
PlotErrorNorm(tgrid,CLsim,PrintFigureTitles)
subplot(3,1,3)
PlotControl(tgrid,CLsim,PrintFigureTitles)

%% State of the controlled PDE

% In plotting and animating the state,
% fill in the homogeneous Dirichlet boundary condition at x=1
spgrid_plot = [spgrid, 1];

figure(3)
colormap jet
Plot1DHeatSurf(CLsim.xesol(1:N,:),spgrid_plot,tgrid,BCtype)

%% Animation of the state of the controlled PDE

figure(4)
% No movie recording
[~,zlims] = Anim1DHeat(CLsim.xesol(1:N,:),spgrid_plot,tgrid,BCtype,0.03,0);

% Movie recording
% [MovAnim,zlims] = Anim1DHeat(CLsim.xesol(1:N,:),spgrid_plot,tgrid,BCtype,0.03,1);

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
