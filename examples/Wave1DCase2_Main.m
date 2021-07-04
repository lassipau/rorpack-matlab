%% Robust output tracking for a 1D wave equation
% with distributed control and observation and possible distributed disturbance input.
% Dirichlet boundary conditions at $x_i=0$ and $x_i=1$.
% Since the control and observation are distributed,
% the system is only strongly (and polynomially) stabilizable.
% Because of this, the low-gain controller cannot be used,
% and the observer-based controller designs do not guarantee closed-loop stability.
% However, since the system is passive,
% the Passive Robust Controller can be used in achieving robust output regulation.

addpath(genpath('../RORPack/'))

% Parameters for this example.
N = 40;
% Construct the system and define the initial state. 
% Input and disturbance input profiles
Bfun = @(x) 10 * (1 - x);
Bdfun = @(x) 5 * x .* (1 - x);
% ('w0' = initial profile, 'wd0' = initial velocity)
w0fun = @(x) zeros(size(x));
% w0fun = @(x) 1 + cos(3 * pi * x) + cos(6*x);
% w0fun = @(x) x .* (x - 1) .* (2 - 5 * x);
wd0fun = @(x) zeros(size(x));

[Sys, x0, phin] = ConstrWave1DCase2(N, Bfun, Bdfun, w0fun, wd0fun);


% Define the reference and disturbance signals:
% Note that the system has an invariant zero at s=0, and therefore the
% regulation of constant signals is impossible

% Case 1:
% yref = @(t) sin(pi*t) + 0.25*cos(2*pi*t);
% wdist = @(t) sin(6*pi*t);
% freqsReal = [pi 2*pi 6*pi];

% Case 2:
yref = @(t) sin(pi*t) + 0.25*cos(2*pi*t);
wdist = @(t) zeros(size(t));
freqsReal = [pi 2*pi];

% Alternative: Tracking of general 2-periodic reference signals
% Tracking of an arbitrary (nonsmooth) 2-periodic signal (no disturbance input)
% wdist = @(t) zeros(size(t));
 
% Begin by defining the function on a single period 0<t<2
% A nonsmooth triangle signal
% yref1per = @(t) (2.*t-1).*(t>=0).*(t<=1)+(3-2.*t).*(t>1).*(t<2);

% Semicircles
% yref1per = @(t) sqrt(1-(t-1).^2);

% Alternating semicircles
% yref1per = @(t) sqrt(abs(1/4-(t-1/2).^2)).*(t>=0).*(t<1)-sqrt(abs(1/4-(t-3/2).^2)).*(t>=1).*(t<2);

% Bump and constant
% yref1per = @(t) sqrt(abs(1/4-(t-1).^2)).*(t>=1/2).*(t<3/2);
% yref1per = @(t) sqrt(abs(1/4-(t-1/2).^2)).*(t>=0).*(t<1);

% The constant part of the signal cannot be tracked due to the second order
% zero of the plant at zero.
% We therefore normalize yref(t) to have average zero on [0,2]
% yr_ave = integral(yref1per,0,2);
% yref = @(t) yref1per(mod(t,2)) - yr_ave/2;

% Include frequencies pi*k that are required in tracking 2-periodic signals
% freqsReal = pi*(1:9);


% Check that the system, the reference and disturbance signals, and the
% frequencies are defined in a consistent way.
Sys = SysConsistent(Sys,yref,wdist,freqsReal);

%% Construct the controller

% A Passive Robust Controller
%
% Since the system is unstable, requires prestabilization with a
% controller feedthrough term 'Dc'. The passivity of the system implies
% that we can use negative output feedback (with any negative definite Dc),
% but the resulting closed-loop system is not exponentially stable, but 
% only strongly (and sometimes polynomially) stable.

dimY = size(Sys.C,1);

% Negative feedback gain for output stabilization
kappa_S = -1;
Dc = kappa_S*eye(dimY);

epsgain = [0.01, 0.3];
% epsgain = 0.3;

[ContrSys, epsgain] = PassiveRC(freqsReal, dimY, epsgain, Sys, Dc);
epsgain


%% Closed-loop simulation and visualization of the results

% Construct the closed-loop system
CLSys = ConstrCLSys(Sys,ContrSys);

stabmarg = CLStabMargin(CLSys)


% Simulate the closed-loop system. 
% The initial state z0 of the controller is chosen to be zero by default.
z0 = zeros(size(ContrSys.G1,1),1);
xe0 = [x0;z0];

Tend = 24;
% Tend = 14;
tgrid = linspace(0, Tend, 601);

CLsim = SimCLSys(CLSys,xe0,yref,wdist,tgrid,[]);


% Choose whther or not to print titles of the figures
PrintFigureTitles = true;
spgrid = linspace(0,1,N);

figure(1)
subplot(3,1,1)
PlotOutput(tgrid,yref,CLsim,PrintFigureTitles)
subplot(3,1,2)
PlotErrorNorm(tgrid,CLsim,PrintFigureTitles)
subplot(3,1,3)
PlotControl(tgrid,CLsim,PrintFigureTitles)

%% Plot the state of the controlled system
figure(3)
colormap jet
surf_plotskip = 2;
Plot1DWaveSurf(CLsim.xesol(1:2*N,1:surf_plotskip:end),phin,spgrid,tgrid(1:surf_plotskip:end))
% Plot1DWaveSurf(CLsim.xesol(1:2*N,:)-CLsim.xesol(dimX+dimZ+(1:(2*N)),:),phin,spgrid,tgrid,[-9 9])
% Plot1DWaveSurf(CLsim.xesol(dimX+dimZ+(1:(2*N)),:),phin,spgrid,tgrid,[-4 4])
set(gca,'ztick',-8:4:8);


%% Plot the reference signal
figure(4)
tt = linspace(0,16,500);
plot(tt,yref(tt),'Color',1.1*[0 0.447 0.741],'Linewidth',3);
set(gca,'xgrid','on','ygrid','on','tickdir','out','box','off')

%% Animation
figure(5)
colormap jet
% No movie recording
[~,zlims] = Anim1DWaveSpectral(CLsim.xesol(1:2*N,:),phin,spgrid,tgrid,0.03,0);
