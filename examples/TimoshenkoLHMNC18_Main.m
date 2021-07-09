%% Robust output tracking of a passive Timoshenko beam model with distributed inputs and outputs. 
% Originally written for the conference presentation at Lagrangian and 
% Hamiltonian Methods in Nonlinear Control (LHMNC) in Valparaiso, Chile in
% 2018. The simulation is associated to the conference paper by Paunonen,
% Le Gorrec, and Ramirez at LHMNC 2018 (the conference paper does not
% include the simulation).
% The controller is a "simple passive controller" studied in the conference
% paper (also see Rebarber-Weiss 2003).
% The system has collocated distributed inputs and outputs, and in the 
% simulation it is approximated using a Finite Differences scheme.


% addpath(genpath('../RORPack/'))


N = 50; 

% Initial state of the plant (the state has four components, x_1, ... x_4)
% The initial state is based on the initial data for w(xi,t) and phi(xi,t),
% namely w0fun (initial deflection), wd0fun (initial velocity), 
% phi0fun (initial angular displacement), phid0fun (initial angular
% velocity)

w0fun = @(x) cos(pi*x);
% wd0fun = @(x) zeros(size(x));
wd0fun = @(x) 5*sin(pi*x);
phi0fun = @(x) zeros(size(x));
phid0fun = @(x) zeros(size(x));


[x0,spgrid,Sys] = ConstrTimoshenkoLHMNC18(w0fun,wd0fun,phi0fun,phid0fun,N);


% Define the reference and disturbance signals
% NOTE: The system has a transmission zero at s=0, and thus the tracking 
% and rejection of signal parts with this frequency is not possible! (i.e.,
% freqsReal should not contain 0).

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


freqsReal = [1, 2];


% Check that the system, the reference and disturbance signals, and the
% frequencies are defined in a consistent way.
Sys = SysConsistent(Sys,yref,wdist,freqsReal);

%% Controller construction

% Simple passive robust controller, used in the original Timoshenko beam
% example in the LHMNC 2018 conference paper (simulation not included in
% the paper).

dimY = size(Sys.C,1);
epsgain = [3,7];
% epsgain = 13;
[ContrSys,epsgain] = PassiveRC(freqsReal,dimY,epsgain,Sys);
epsgain

% % Alternative controller:
% % An observer-based robust controller
% % Stabilizing state feedback and output injection operators K and L
% % These are chosen based on collocated design. 
% K21 = -.5*Sys.B';
% % PlotEigs(full(Sys.A+Sys.B*K),[-1 .1 NaN NaN])
% 
% L = -.5*Sys.C';
% % PlotEigs(full(Sys.A+L*Sys.C))
% 
% IMstabtype = 'poleplacement';
% % IMstabtype = 'LQR';
% IMstabmarg = 0.5;
% 
% ContrSys = ObserverBasedRC(freqsReal,Sys,K21,L,IMstabtype,IMstabmarg);

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
PlotEigs(full(CLSys.Ae))

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

% In order to visualize the behaviour of the beam, we need to compute the
% deflection profile w(xi,t) based on the state variable x_2(t) =
% \rho*\dot{w}(\xi,t) by numerical integration. 
% We use a denser grid for 't' for the numerical integration
tt_int = linspace(0,Tend,601);
xe_int = deval(CLsim.solstruct,tt_int);
% In this example rho=1
rho = 1;
profile_int = 1/rho*(w0fun(spgrid)*ones(1,length(tt_int))+[zeros(1,length(tt_int));cumtrapz(tt_int(2)-tt_int(1),xe_int((N+1):(2*N),:),2)]);
% Interpolate the data to the plotting grid
profile = interp1(tt_int.',profile_int.',tgrid(:)).';


figure(3)
plotskip = 2;
PlotLHMNCSurf(profile(:,1:plotskip:end),spgrid,tgrid(:,1:plotskip:end))


%% The reference signal

figure(4)
tt = linspace(0,16,500);
plot(tt,yref(tt),'Color',1.1*[0 0.447 0.741],'Linewidth',3);
title('Reference signal $y_{ref}$','Interpreter','latex','Fontsize',16)
set(gca,'xgrid','on','ygrid','on','tickdir','out','box','off')
