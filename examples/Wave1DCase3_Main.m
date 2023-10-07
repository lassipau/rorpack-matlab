%% A Wave equation in 1D with Neumann boundary control and pointwise observation, 
% The system has viscous (spatially varying) damping, Neumann boundary
% input and input disturbance, and pointwise measurement (either deflection
% or velocity at a given point).
% Simulation based on modal approximation.
% 
% In this example, the system is exponentially and has noncollocated I/O. 

% addpath(genpath('../RORPack/'))

N = 60; 


% Define the viscous damping function 'dfun'
% dfun = @(x) zeros(size(x));
dfun = @(x) 20*x.*(1-x);

% Output is the measurement at point xi0, either deflection (p_ind = 1) or
% velocity (p_ind = 2);
xi0 = 0;
p_ind = 1;
% p_ind = 2;



% Initial state of the plant
w0fun = @(x) zeros(size(x));
w0fun = @(x) cos(pi*1.5*x);
wd0fun = @(x) zeros(size(x));
%x0fun = @(x,y) 0.5*(1+cos(pi*(1-x))).*(1-1/4*cos(2*pi*y));
% x0fun = @(x,y) 0.5*(1+cos(pi*(1-x)));
%x0fun = @(x,y) 1/2*x.^2.*(3-2*x)-1;
%x0fun = @(x,y) 1/2*x.^2.*(3-2*x)-1/2;
%x0fun = @(x,y) 1*(1-x).^2.*(3-2*(1-x))-1;
%x0fun = @(x,y) .5*(1-x).^2.*(3-2*(1-x))-.5;
%x0fun = @(x,y) 1/4*(x.^3-1.5*x.^2)-1/4;
%x0fun = @(x,y) .2*x.^2.*(3-2*x)-.5;


[x0,Sys,phin] = ConstrWave1DCase3(dfun,xi0,w0fun,wd0fun,N);

% Select the correct measurement based on 'p_ind'. Either deflection (p_ind = 1) or
% velocity (p_ind = 2);
% Sys.A = Sys.A-.5*Sys.B*(Sys.B)'-.1*(Sys.C)'*Sys.C;
C = Sys.C;
D = Sys.D;
Dd = Sys.Dd;

Sys.C = C(p_ind,:);
Sys.D = D(p_ind,:);
Sys.Dd = Dd(p_ind,:);

% Define the reference and disturbance signals, and list the
% required frequencies in 'freqsReal'

% yref = @(t) sin(2*t)+.1*cos(6*t);
%yref = @(t) sin(2*t)+.2*cos(3*t);
%yref = @(t) ones(size(t));

%wdist = @(t) ones(size(t));
%wdist = @(t) sin(t);
%wdist = @(t) zeros(size(t));

% Triangle signal case
% Begin by defining the function on a single period 0<t<2
% A nonsmooth triangle signal
% yref1per = @(t) (2.*t-1).*(t>=0).*(t<=1)+(3-2.*t).*(t>1).*(t<2);

% Semicircles
% yref1per = @(t) sqrt(1-(t-1).^2);

% Alternating semicircles
% yref1per = @(t) sqrt(abs(1/4-(t-1/2).^2)).*(t>=0).*(t<1)-sqrt(abs(1/4-(t-3/2).^2)).*(t>=1).*(t<2);

% % Bump and constant
% yref1per = @(t) sqrt(abs(1/4-(t-1).^2)).*(t>=1/2).*(t<3/2);
% % yref1per = @(t) sqrt(abs(1/4-(t-1/2).^2)).*(t>=0).*(t<1);
% 
% % The constant part of the signal cannot be tracked due to the second order
% % zero of the plant at zero.
% % We therefore normalize yref(t) to have average zero on [0,2]
% yr_ave = integral(yref1per,0,2);
% yref = @(t) yref1per(mod(t,2)) - yr_ave/2;
% 
% freqsReal = pi*(1:29);


% Case 1:
yref = @(t) .5*sin(2*t)+.5*cos(3*t);
% yref = @(t) sin(pi/2*t)-2*cos(pi/2*t);
% yref = @(t) zeros(size(t));
wdist = @(t) zeros(size(t));
% wdist = @(t) sin(2*pi*t);

% Case 2:
% yref = @(t) ones(size(t));
% wdist = @(t) zeros(size(t));

% Case 3:
% yref = @(t) sin(2*t)+.1*cos(6*t);
% wdist = @(t) sin(t);

freqsReal = [2,3];


% Check that the system, the reference and disturbance signals, and the
% frequencies are defined in a consistent way.
Sys = SysConsistent(Sys,yref,wdist,freqsReal);


%% Construct the controller

% A Low-Gain 'Minimal' Robust Controller

dimX = size(Sys.A,1);
dimY = size(Sys.C,1);
Pappr = @(s) Sys.C*((s*eye(dimX)-Sys.A)\Sys.B)+Sys.D;
Dc = .4; % Prestabilization with negative output feedback
Dc = 0;

% Approximate the transfer function values of the system under the negative
% output feedback. NOTE: If zero frequency is considered, then the transfer
% function at \lambda=0 has to be done using a different approach!
Pvals = cell(1,length(freqsReal));
for ind = 1:length(freqsReal)
    Pval_nominal = Pappr(1i*freqsReal(ind));
    Pvals{ind} = (eye(dimY)-Pval_nominal*Dc)\Pval_nominal;
end

epsgain = [0.01,.6];
% epsgain = .1;
[ContrSys,epsgain] = LowGainRC(freqsReal,Pvals,epsgain,Sys,Dc);
epsgain

% An observer-based robust controller 
% % Step 1: Compute stabilizing operators K21 and L so that A+B*K21 and A+L*C
% % are exponentially stable
% K21 = -lqr(full(Sys.A),Sys.B,eye(size(Sys.A)),1*eye(dimU),zeros(dimX,dimU));
% L = -lqr(full((Sys.A).'),(Sys.C).',eye(size(Sys.A)),0.1*eye(dimY),zeros(dimX,dimY)).';

% K0 = -1/3*[3,2,zeros(1,2*N-2)];
% L0 = -.6*[3;2;zeros(2*N-2,1)];
% L0 = -.3*[0;1;zeros(2*N-2,1)];

% % Prestabilization
% k_m = .3;
% Sys.A = Sys.A-k_m*Sys.B*Sys.Cm;
% 
% kappa_K = .9;
% kappa_L = .8;
% % K_S = -kappa_K*(Sys.B)'+K0;
% % K_S = -kappa_K*Kinf+K0;
% K_S = -kappa_K*Kinf;
% % PlotEigs(Sys.A+Sys.B*K_S,[-1.5, 0.1, NaN, NaN])
%  
% 
% % L = [zeros(dimX,1),-kappa_L*(Sys.Cm)'+L0];
% % L = [L0,-kappa_L*(Sys.Cm)'];
% L = -kappa_L*Linf;
% % L = L0;
% 
% % PlotEigs(Sys.A+L*[Sys.C;Sys.Cm],[-1 0 -10 10])
% % PlotEigs(Sys.A+L*[Sys.C;Sys.Cm])
% % PlotEigs(Sys.A+L*Sys.C,[-1.5, .01, NaN, NaN])
% % PlotEigs(Sys.A+L*Sys.C)
% 
% 
% % [ContrSys,K21] = ObserverBasedRC(freqsReal,Sys,K_S,L,'poleplacement',3);
% [ContrSys,K21] = ObserverBasedRC(freqsReal,Sys,K_S,L,'LQR',3);
% % [ContrSys,K21] = DualObserverBasedRC(freqsReal,Sys,K_S,L,'LQR',3);
% % [ContrSys,K21] = DualObserverBasedRC(freqsReal,Sys,K_S,L,'poleplacement',3);

%% Closed-loop construction and simulation

% Construct the closed-loop system
CLSys = ConstrCLSys(Sys,ContrSys);

% Print an approximate stability margin of the closed-loop system
stabmarg = CLStabMargin(CLSys)


% Define the initial state of the closed-loop system
% (the controller has zero initial state by default).
xe0 = [x0;zeros(size(ContrSys.G1,1),1)];

% Set the simulation length and define the plotting grid
Tend = 20;
tgrid = linspace(0,Tend,401);

% Simulate the closed-loop system
CLsim = SimCLSys(CLSys,xe0,yref,wdist,tgrid,[]);

%% Visualization

% Plot the (approximate) eigenvalues of the closed-loop system
figure(1)
PlotEigs(full(CLSys.Ae),[-2 NaN NaN NaN])

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

% Plots for additional analysis
% q = length(freqsReal);
% if freqsReal(1)==0,dimZ = 2*q-1; else dimZ=2*q;end
% dimX = size(Sys.A,1);
% dimU = size(ContrSys.K,1);
% K1full = [ContrSys.K(:,1:dimZ), zeros(dimU,dimX)];
% K2full = [zeros(dimU,dimZ) ContrSys.K(:,(dimZ+1):end)];
% figure(6)
% plot(tgrid,[CLsim.control;K1full*CLsim.xesol((2*N+1):end,:);K2full*CLsim.xesol((2*N+1):end,:)],'Linewidth',2);
% obserror = sum((CLsim.xesol(1:2:(2*N),:)-CLsim.xesol(dimX+dimZ+(1:2:(2*N)),:)).^2,1);
% figure(7)
% plot(tgrid,obserror,'Linewidth',2);
% title('Observer error in the controller.')

%% State of the controlled PDE

figure(3)
spgrid = linspace(0,1,41);
% To use Plot1DWaveSurf, we need to reorganise the indices of the state
% variables
index_perm = [1:N;(N+1):(2*N)];
index_perm = index_perm(:);
AnimData = CLsim.xesol(index_perm,:);
Plot1DWaveSurf(AnimData,phin,spgrid,tgrid)
% Plot1DWaveSurf(CLsim.xesol(1:2*N,:)-CLsim.xesol(dimX+dimZ+(1:(2*N)),:),phin,spgrid,tgrid,[-9 9])
% Plot1DWaveSurf(CLsim.xesol(dimX+dimZ+(1:(2*N)),:),phin,spgrid,tgrid,[-4 4])
% set(gca,'ztick',-8:4:8);
% colormap jet

%% Animation

figure(4)
% No movie recording
[~,zlims] = Anim1DWaveSpectral(AnimData,phin,spgrid,tgrid,0.01,0);

% Tpause = 0.05;
% record = 0;
% [MovAnim,zlims] = Anim1DWaveSpectral(CLsim.xesol(1:2*N,:),phin,spgrid,tgrid,Tpause,record)

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