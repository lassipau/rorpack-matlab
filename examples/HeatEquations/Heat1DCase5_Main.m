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

N = 100; 

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

[x0,Sys,spgrid,BCtype] = ConstrHeat1DCase5(cfun,x0fun,N);

% Define the reference and disturbance signals (the system has two outputs
% and two disturbance inputs)
% Case 1:
yref = @(t) [sin(2*t);2*cos(3*t)];
wdist = @(t) [sin(6*t);-ones(size(t))];
% wdist = @(t) zeros(2,size(t));

% Case 2:
% yref = @(t) [sin(2*t)+.3*cos(6*t);ones(size(t))];
% wdist = @(t) [sin(1*t);ones(size(t))+.3*sin(3*t+0.3)];


freqsReal = [1, 2, 3, 6];

% Check the consistency of the system definition
Sys = SysConsistent(Sys,yref,wdist,freqsReal);


%% Construct the controller

% % A Low-Gain 'Minimal' Robust Controller
%
% Pappr = @(s) Sys.C*((s*eye(size(Sys.A,1))-Sys.A)\Sys.B)+Sys.D;
% Pvals = cell(1,length(freqs));
% for ind = 1:length(freqs)
%   Pvals{ind} = Pappr(freqs(ind));
% end
% 
% epsgain = [0.01,6];
% % epsgain = .1;
% [ContrSys,epsgain] = LowGainRC(freqs,Pvals,epsgain,Sys);
% epsgain

% A Passive Robust Controller
% Since the system is unstable, requires prestabilization with a
% controller feedthrough term.

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
% % K = -10*Sys.B';
% K = -lqr(Sys.A,Sys.B,0.1*eye(N),10*eye(2));
% % K = zeros(2,N);
% % PlotEigs(full(Sys.A+Sys.B*K),[NaN .1 -.3 .3])
% 
% % L = -10*Sys.C';
% L = -lqr(Sys.A',Sys.C',10*eye(N),eye(2))';
% % PlotEigs(full(Sys.A+L*Sys.C),[-20 1 -.3 .3])
% 
% ContrSys = ObserverBasedRC(freqsReal,Sys,K,L,'LQR',1);
% % ContrSys = ObserverBasedRC(freqsReal,Sys,K,L,'poleplacement',1);
% % ContrSys = DualObserverBasedRC(freqsReal,Sys,K,L,'LQR',1);
% % ContrSys = DualObserverBasedRC(freqsReal,Sys,K,L,'poleplacement',1);
 

% % A reduced order observer-based robust controller
% % DISCLAIMER: This is not directly supported by the theoretical results in
% % Paunonen-Phan 2020 due to the control at the boundaries!
% %
% % The construction of the controller uses a Galerkin approximation
% % of the heat system: The Galerkin arpproximation used in the controller
% % design is a lower dimensional numerical approximation of the PDE model.
% Nlow = 50;
% [~,Sys_Nlow,~,~] = Constr1DHeatCase5(cfun,x0fun,Nlow);
% 
% % Store the Galerkin approximation in "SysApprox".
% SysApprox.AN = Sys_Nlow.A;
% SysApprox.BN = Sys_Nlow.B;
% SysApprox.CN = Sys_Nlow.C;
% SysApprox.D = Sys_Nlow.D;
% alpha1 = 1;
% alpha2 = 0.5;
% Q0 = eye(IMdim(freqs,size(SysApprox.CN,1))); % Size = dimension of the IM 
% Q1 = eye(size(SysApprox.AN,1)); % Size = dim(V_N)
% Q2 = eye(size(SysApprox.AN,1)); % Size = dim(V_N)
% R1 = eye(size(SysApprox.CN,1)); % Size = dim(Y)
% R2 = eye(size(SysApprox.BN,2)); % Size = dim(U)
% ROMorder = 3;
% 
% ContrSys = ObserverBasedROMRC(freqsReal,SysApprox,alpha1,alpha2,R1,R2,Q0,Q1,Q2,ROMorder);


%% Closed-loop simulation
CLSys = ConstrCLSys(Sys,ContrSys);

stabmarg = CLStabMargin(CLSys)

figure(1)
PlotEigs(CLSys.Ae,[-30 .3 -6 6])


xe0 = [x0;zeros(size(ContrSys.G1,1),1)];

Tend = 8;
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
