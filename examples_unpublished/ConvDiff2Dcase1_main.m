%% Heat equation on a square with Neumann boundary conditions and 
% Neumann boundary control and Dirichlet boundary observation 
% Approximation with a Finite differences scheme 

% Unstable system, stabilization by stabilizing the only unstable
% eigenvalue =0

% addpath(genpath('../RORPack/'))




h_param = 0.08;
% h_param = 0.1;
refine_times = 1;

% Initial state of the plant (function of two variables)
%x0fun = @(x,y) zeros(size(x));
x0fun = @(x,y) 0.5*(1+cos(pi*(1-x))).*(1-1/4*cos(2*pi*y));
% x0fun = @(x,y) 0.5*(1+cos(pi*(1-x)));
%x0fun = @(x,y) 1/2*x.^2.*(3-2*x)-1;
%x0fun = @(x,y) 1/2*x.^2.*(3-2*x)-1/2;
%x0fun = @(x,y) 1*(1-x).^2.*(3-2*(1-x))-1;
%x0fun = @(x,y) .5*(1-x).^2.*(3-2*(1-x))-.5;
%x0fun = @(x,y) 1/4*(x.^3-1.5*x.^2)-1/4;
%x0fun = @(x,y) .2*x.^2.*(3-2*x)-.5;

% % The spatially varying thermal diffusivity of the material
% % cfun = @(t) ones(size(t));
% % cfun = @(t) 1+t;
% cfun = @(t) 1-2*t.*(1-2*t);
% % cfun = @(t) 1+0.5*cos(5/2*pi*t);
% % cfun = @(t) 0.3-0.6*t.*(1-t);

% [x0,spgrid,Sys] = ConstrHeat2Dtest1(1,x0fun,N);

testcase = 4;

[x0,Sys] = Constr2DConvDiff_case1(h_param,refine_times,testcase,x0fun);


%%

%yref = @(t) sin(2*t)+.1*cos(6*t);
%yref = @(t) sin(2*t)+.2*cos(3*t);
%yref = @(t) ones(size(t));

%wdist = @(t) ones(size(t));
%wdist = @(t) sin(t);
%wdist = @(t) zeros(size(t));

% Case 1:
yref = @(t) [sin(2*t);-2*ones(size(t))];%+.2*cos(3*t);
%wdist = @(t) zeros(size(t));
wdist = @(t) sin(2*t);

% Case 2:
% yref = @(t) ones(size(t));
% wdist = @(t) ones(size(t));

% Case 3:
% yref = @(t) sin(2*t)+.1*cos(6*t);
% wdist = @(t) sin(t);

freqsReal = [0, 1, 2, 3, 6];

% Check the consistency of the system definition
Sys = SysConsistent(Sys,yref,wdist,freqsReal);

%%

% A reduced order observer-based robust controller
%
% The construction of the controller uses a Galerkin approximation
% of the heat system: The Galerkin arpproximation used in the controller
% design is a lower dimensional numerical approximation of the PDE model.
% This can be achieved using a larger mesh parameter 'h_param', and
% possibly smaller number of refinement interations ('refine_times')


h_param = 0.085;
% h_param = 0.1;
refine_times = 1;
[~,Sys_Nlow] = Constr2DConvDiff_case1(h_param,refine_times,testcase,x0fun);


% Store the numerical approximation in "SysApprox".
SysApprox.AN = Sys_Nlow.A;
SysApprox.BN = Sys_Nlow.B;
SysApprox.CN = Sys_Nlow.C;
SysApprox.D = Sys_Nlow.D;
alpha1 = 0.4;
alpha2 = 0.4;
Q0 = eye(IMdim(freqsReal,size(SysApprox.CN,1))); % Size = dimension of the IM 
Q1 = eye(size(SysApprox.AN,1)); % Size = dim(V_N)
Q2 = eye(size(SysApprox.AN,1)); % Size = dim(V_N)
R1 = 0.3*eye(size(SysApprox.CN,1)); % Size = dim(Y)
R2 = 0.3*eye(size(SysApprox.BN,2)); % Size = dim(U)
ROMorder = 50;
% ROMorder = NaN; % No model reduction for the controller

ContrSys = ObserverBasedROMRC(freqsReal,SysApprox,alpha1,alpha2,R1,R2,Q0,Q1,Q2,ROMorder);


CLSys = ConstrCLSys(Sys,ContrSys);

% Print the stability margin, stop the simulation run if the closed-loop 
% system is unstable
stabmarg = CLStabMargin(CLSys,true)

figure(1)
PlotEigs(CLSys.Ae); %,[-1 .3 -4 4]);

% PlotEigs(CLSys.Ae,[-2 .3 -6 6]);
%%


xe0 = [x0;zeros(size(ContrSys.G1,1),1)];

Tend = 30;
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

% ANIMATION NOT YET IMPLEMENTED

% figure(3)
% colormap jet
% % No movie recording
% [~,zlims] = AnimHeat2Dtest1(CLsim,spgrid,tgrid,0.03,0);
% 
% % Movie recording
% % [MovAnim,zlims] = AnimHeat2Dtest1(CLsim,spgrid,tgrid,0,1);
% 
% %movie(MovAnim)

%%

% figure(4)
% colormap jet
% %PlotHeat2DSurf(x0,spgrid,[-1.4,1.4])
% PlotHeat2DSurf(x0,spgrid,zlims)
% 
% figure(5)
% tt = linspace(0,16,500);
% plot(tt,yref(tt),'Color',1.1*[0 0.447 0.741],'Linewidth',3);
% title('Reference signal $y_{ref}$','Interpreter','latex','Fontsize',16)
% set(gca,'xgrid','on','ygrid','on','tickdir','out','box','off')


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
