%% Heat equation on the interval [0,1] with 
% Neumann boundary control and Dirichlet boundary observation 
% Approximation with a Finite differences scheme 

% Case 1: Neumann boundary control at x=0, regulated output y(t) and a 
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

[x0,Sys,spgrid,BCtype] = Constr1DHeatCase2(1,x0fun,N);


% Model = ss(Sys.A,Sys.B,Sys.C,Sys.D);
% tt=linspace(0,4);
% [output,t,xx]=lsim(Model,ones(size(tt)),tt);
% % plot(tt,(xx(:,2)-xx(:,1))*(N-1))
% surf(linspace(0,1,50),tt,xx)

%%

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
freqs = [0 1 2 3 6];

% Sys.A = Sys.A+2*pi^2*Sys.B*Sys.Cm;
% PlotEigs(full(Sys.A),[-20 1 -.3 .3])

% eig(full(Sys.A))

% A Low-Gain 'Minimal' Robust Controller

% dimX = size(Sys.A,1);
% Pappr = @(s) Sys.C*((s*eye(dimX)-Sys.A)\Sys.B)+Sys.D;
% 
% Pvals = cell(1,length(freqs));
% for ind = 1:length(freqs)
%   Pvals{ind} = Pappr(freqs(ind));
% end
% 
% epsgainrange = [0.01,3];
% % epsgain = .1;
% [ContrSys,epsgain] = ConstrContrLG(freqs,Pvals,epsgainrange,Sys);
% epsgain

% An observer-based robust controller
% Stabilizing state feedback and output injection operators K and L
% These are chosen based on collocated design. Only the single unstable
% eigenvalue at s=0 needs to be stabilized
K = 7*[1, zeros(1,N-1)];
% K = zeros(1,N);
% PlotEigs(full(Sys.A+Sys.B*K),[NaN .1 -.3 .3])

L = -7*[2*N;zeros(N-1,1)];
% PlotEigs(full(Sys.A+L*Sys.C),[-20 1 -.3 .3])

ContrSys = ConstrContrObsBased(freqs,Sys,K,L,'LQR',0.5);
% ContrSys = ConstrContrObsBased(freqs,Sys,K,L,'poleplacement',0.5);
% ContrSys = ConstrContrDualObsBased(freqs,Sys,K,L,'LQR',0.5);
% ContrSys = ConstrContrDualObsBased(freqs,Sys,K,L,'poleplacement',0.5);

%% Closed-loop simulation
CLSys = ConstrCLSys(Sys,ContrSys);

stabmarg = CLStabMargin(CLSys)

figure(1)
PlotEigs(CLSys.Ae,[-20 .3 -6 6])
%%

xe0 = [x0;zeros(size(ContrSys.G1,1),1)];

Tend = 14;
tgrid = linspace(0,Tend,300);


CLsim = SimCLSys(CLSys,xe0,yref,wdist,tgrid,[]);

% Choose whther or not to print titles of the figures
PrintFigureTitles = true;

figure(2)
subplot(3,1,1)
plotOutput(tgrid,yref,CLsim,PrintFigureTitles)
subplot(3,1,2)
plotErrorNorm(tgrid,CLsim,PrintFigureTitles)
subplot(3,1,3)
plotControl(tgrid,CLsim,ContrSys,N,PrintFigureTitles)

%%


figure(3)
colormap jet
Plot1DHeatSurf(CLsim.xesol(1:N,:),spgrid,tgrid,BCtype)

%%
% figure(4)
% No movie recording
% [~,zlims] = Anim1DHeat(CLsim.xesol(1:N,:),spgrid,tgrid,BCtype,0.03,0);

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
