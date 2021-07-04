%% Heat equation on the interval [0,1] with 
% Neumann boundary control and Dirichlet boundary observation 
% Approximation with a Finite differences scheme 

% Neumann boundary disturbance at x=0,
% two distributed controls and
% two distributed measurement regulated outputs. The controls act 
% on the intervals 'IB1' and 'IB2' (Default 'IB1' = [0.3,0.4] and 
% 'IB2' = [0.6,0.7]) and the measurements are the average 
% temperatures on the intervals 'IC1'  and 'IC2' (Default 
% 'IC1' = [0.1,0.2] and 'IC2' = [0.8,0.9]).

addpath(genpath('../RORPack/'))

% Parameters for this example.
N = 100; 

% Initial state of the plant
%x0fun = @(x) zeros(size(x));
x0fun = @(x) 1*(1+cos(pi*(1-x)));
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
% cfun = @(t) 1+0.5*cos(5/2*pi*t);
cfun = @(t) 0.3-0.6*t.*(1-t);

IB1 = [.3, .4];
IB2 = [.6, .7];
IC1 = [.1, .2];
IC2 = [.8, .9];

[x0,Sys,spgrid,BCtype] = ConstrHeat1DCase3(cfun,x0fun,N,IB1,IB2,IC1,IC2);


% Model = ss(Sys.A,Sys.B,Sys.C,Sys.D);
% tt=linspace(0,4);
% [output,t,xx]=lsim(Model,ones(2,length(tt)),tt,ones(N,1));
% h = spgrid(2)-spgrid(1);
% % plot(tt,output,'Linewidth',2)
% plot(tt,sum(h*xx.^2,2),'Linewidth',2)
% % plot(tt,(xx(:,2)-xx(:,1))*(N-1))

% Define the reference and disturbance signals, and list the
% required frequencies in 'freqsReal'

%yref = @(t) sin(2*t)+.1*cos(6*t);
%yref = @(t) sin(2*t)+.2*cos(3*t);
%yref = @(t) ones(size(t));

%wdist = @(t) ones(size(t));
%wdist = @(t) sin(t);
%wdist = @(t) zeros(size(t));

% Case 1:
% yref = @(t) [sin(2*t);2*cos(3*t)];
% wdist = @(t) zeros(size(t));
% wdist = @(t) sin(6*t);
% freqsReal = [2 3 6];

% Case 2:
% yref = @(t) ones(2, size(t,1));
% wdist = @(t) zeros(size(t));
yref = @(t) [2*ones(size(t));2*cos(3*t)+sin(2*t)];
wdist = @(t) 6*cos(t);
freqsReal = [1 2 3 6];

% Sys.A = Sys.A-Sys.B*Sys.B';
% PlotEigs(full(Sys.A),[-20 1 -.3 .3])

% eig(full(Sys.A))

% Check the consistency of the system definition
Sys = SysConsistent(Sys,yref,wdist,freqsReal);


%% Construct the controller

% A Low-Gain 'Minimal' Robust Controller

% Pappr = @(s) Sys.C*((s*eye(size(Sys.A,1))-Sys.A)\Sys.B)+Sys.D;
% 
% Pvals = cell(1,length(freqsReal));
% for ind = 1:length(freqsReal)
%   Pvals{ind} = Pappr(freqsReal(ind));
% end
% 
% epsgainrange = [0.01,6];
% % epsgain = .1;
% [ContrSys,epsgain] = LowGainRC(freqsReal,Pvals,epsgainrange,Sys);
% epsgain

% An observer-based robust controller
% Stabilizing state feedback and output injection operators K and L
% These are chosen based on collocated design. Only the single unstable
% eigenvalue at s=0 needs to be stabilized
K = -Sys.B';
% K = zeros(2,N);
% PlotEigs(full(Sys.A+Sys.B*K),[NaN .1 -.3 .3])

%
L = -10*Sys.C';
%  PlotEigs(full(Sys.A+L*Sys.C),[-20 1 -.3 .3])

% ContrSys = ObserverBasedRC(freqsReal,Sys,K,L,'LQR',1);
% ContrSys = ObserverBasedRC(freqsReal,Sys,K,L,'poleplacement',1);
% ContrSys = DualObserverBasedRC(freqsReal,Sys,K,L,'LQR',1);
% ContrSys = DualObserverBasedRC(freqsReal,Sys,K,L,'poleplacement',1);

% A reduced order observer-based robust controller
%
% The construction of the controller uses a Galerkin approximation
% of the heat system:
% The Galerkin arpproximation used in the controller
% design is a lower dimensional numerical approximation
% of the PDE model.
Nlow = 50;
[~,Sys_Nlow,~,~] = ConstrHeat1DCase3(cfun,x0fun,Nlow,IB1,IB2,IC1,IC2);

% Store the numerical approximation in "SysApprox".
SysApprox.AN = Sys_Nlow.A;
SysApprox.BN = Sys_Nlow.B;
SysApprox.CN = Sys_Nlow.C;
SysApprox.D = Sys_Nlow.D;
alpha1 = 1.5;
alpha2 = 1;
Q0 = eye(IMdim(freqsReal,size(SysApprox.CN,1))); % Size = dimension of the IM 
Q1 = eye(size(SysApprox.AN,1)); % Size = dim(V_N)
Q2 = eye(size(SysApprox.AN,1)); % Size = dim(V_N)
R1 = eye(size(SysApprox.CN,1)); % Size = dim(Y)
R2 = eye(size(SysApprox.BN,2)); % Size = dim(U)
ROMorder = 3;

ContrSys = ObserverBasedROMRC(...
    freqsReal,SysApprox,alpha1,alpha2,R1,R2,Q0,Q1,Q2,ROMorder);

Sys = SysConsistent(Sys,yref,wdist,freqsReal);
%% Closed-loop construction and simulation

% Construct the closed-loop system
CLSys = ConstrCLSys(Sys,ContrSys);

stabmarg = CLStabMargin(CLSys)

xe0 = [x0;zeros(size(ContrSys.G1,1),1)];

Tend = 8;
tgrid = linspace(0,Tend,300);

CLsim = SimCLSys(CLSys,xe0,yref,wdist,tgrid,[]);

%% Visualization

figure(1)
PlotEigs(CLSys.Ae,[-30 .3 -6 6])

% Choose whether or not to print titles of the figures
PrintFigureTitles = true;

figure(2)
subplot(3,1,1)
PlotOutput(tgrid,yref,CLsim,PrintFigureTitles)
subplot(3,1,2)
PlotErrorNorm(tgrid,CLsim,PrintFigureTitles)
subplot(3,1,3)
PlotControl(tgrid,CLsim,PrintFigureTitles)


% In plotting and animating the state,
% fill in the homogeneous Dirichlet boundary condition at x=1
spgrid_plot = [spgrid, 1];

figure(3)
colormap jet
surf_t_plotskip = 1;
surf_n_plotskip = 2;
Plot1DHeatSurf(CLsim.xesol(1:surf_n_plotskip:N,1:surf_t_plotskip:end),spgrid_plot(1:surf_n_plotskip:end),tgrid(1:surf_t_plotskip:end),BCtype)

%%
figure(4)
% No movie recording
[~,zlims] = Anim1DHeat(CLsim.xesol(1:N,:),spgrid_plot,tgrid,BCtype,0.03,0);

% Movie recording
% [MovAnim,zlims] = Anim1DHeat(CLsim.xesol(1:N,:),spgrid,tgrid,BCtype,0.03,1);

%movie(MovAnim)

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
