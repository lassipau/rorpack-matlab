%% Heat equation on the interval [0,1] with 
% Neumann boundary control and Dirichlet boundary observation 
% Approximation with a Finite differences scheme 

% Case 1: Neumann boundary control at x=0, regulated output y(t) and a 
% Neumann boundary disturbance at x=1
% Unstable system, stabilization by stabilizing the only unstable
% eigenvalue =0

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

IB1 = [.3, .4];
IB2 = [.6, .7];
IC1 = [.1, .2];
IC2 = [.8, .9];

[x0,Sys,spgrid,BCtype] = Constr1DHeatCase5(x0fun,N);



% Model = ss(Sys.A,Sys.B,Sys.C,Sys.D);
% tt=linspace(0,4);
% [output,t,xx]=lsim(Model,ones(2,length(tt)),tt,ones(N,1));
% h = spgrid(2)-spgrid(1);
% % plot(tt,output,'Linewidth',2)
% plot(tt,sum(h*xx.^2,2),'Linewidth',2)
% % plot(tt,(xx(:,2)-xx(:,1))*(N-1))



%

%yref = @(t) sin(2*t)+.1*cos(6*t);
%yref = @(t) sin(2*t)+.2*cos(3*t);
%yref = @(t) ones(size(t));

%wdist = @(t) ones(size(t));
%wdist = @(t) sin(t);
%wdist = @(t) zeros(size(t));

% Case 1:
yref = @(t) [sin(2*t);2*cos(3*t)];
wdist = @(t) zeros(size(t));
% wdist = @(t) sin(6*t);

% Case 2:
% yref = @(t) ones(size(t));
% wdist = @(t) ones(size(t));

% Case 3:
% yref = @(t) sin(2*t)+.1*cos(6*t);
% wdist = @(t) sin(t);


% freqs = [-3i -2i -1i 0 1i 2i 3i];
freqsReal = [1 2 3 6];
freqs = 1i*freqsReal;


% Sys.A = Sys.A-Sys.B*Sys.B';
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
% epsgainrange = [0.01,6];
% % epsgain = .1;
% [ContrSys,epsgain] = ConstrContrLGReal(freqsReal,Pvals,epsgainrange,Sys);
% epsgain

% An observer-based robust controller
% Stabilizing state feedback and output injection operators K and L
% These are chosen based on collocated design. Only the single unstable
% eigenvalue at s=0 needs to be stabilized
% K = -10*Sys.B';
K = -lqr(Sys.A,Sys.B,0.1*eye(N),10*eye(2));
% K = zeros(2,N);
% PlotEigs(full(Sys.A+Sys.B*K),[NaN .1 -.3 .3])


%
L = -10*Sys.C';
L = -lqr(Sys.A',Sys.C',10*eye(N),eye(2))';
PlotEigs(full(Sys.A+L*Sys.C),[-20 1 -.3 .3])

ContrSys = ConstrContrObsBasedReal(freqsReal,Sys,K,L,'LQR',1);
% ContrSys = ConstrContrObsBasedReal(freqsReal,Sys,K,L,'poleplacement');


%% Closed-loop simulation
CLSys = ConstrCLSys(Sys,ContrSys);

stabmarg = CLStabMargin(CLSys)

PlotEigs(CLSys.Ae,[-30 .3 -6 6]);
%%


xe0 = [x0;zeros(size(ContrSys.G1,1),1)];

Tend = 8;
tgrid = linspace(0,Tend,300);



CLsim = SimCLSys(CLSys,xe0,yref,wdist,tgrid,[]);

figure(1)
subplot(3,1,1)
hold off
cla
hold on
plot(tgrid,yref(tgrid),'Color',1.1*[0 0.447 0.741],'Linewidth',2);
plot(tgrid,CLsim.output,'Color', [0.85 0.325 0.098],'Linewidth',2);
title('Output $y(t)$ (red) and the reference $y_{ref}(t)$ (blue)','Interpreter','latex','Fontsize',16)
set(gca,'xgrid','off','tickdir','out','box','off')
subplot(3,1,2)
plot(tgrid,CLsim.error,'Linewidth',2);
set(gca,'xgrid','on','ygrid','on','tickdir','out','box','off')
title('Regulation error $y(t)-y_{ref}(t)$','Interpreter','latex','Fontsize',16)
set(gcf,'color',1/255*[252 247 255])
subplot(3,1,2)
plot(tgrid,CLsim.error,'Linewidth',2);
set(gca,'xgrid','on','ygrid','on','tickdir','out','box','off')
title('Regulation error $y(t)-y_{ref}(t)$','Interpreter','latex','Fontsize',16)
set(gcf,'color',1/255*[252 247 255])
subplot(3,1,3)
% Plot the control input
plot(tgrid,[zeros(size(ContrSys.K,1),N),ContrSys.K]*CLsim.xesol,'Linewidth',2);
set(gca,'xgrid','on','ygrid','on','tickdir','out','box','off')
title('Control input $u(t)$','Interpreter','latex','Fontsize',16)
set(gcf,'color',1/255*[252 247 255])


%%


figure(2)
colormap jet
Plot1DHeatSurf(CLsim.xesol(1:N,:),spgrid,tgrid,BCtype)

%%
figure(3)
% No movie recording
[~,zlims] = Anim1DHeat(CLsim.xesol(1:N,:),spgrid,tgrid,BCtype,0.03,0);

% Movie recording
%  [MovAnim,zlims] = Anim1DHeat(CLsim.xesol(1:N,:),spgrid,tgrid,BCtype,0.03,1);

% movie(MovAnim)

%%


figure(4)
tt = linspace(0,16,500)
plot(tt,yref(tt),'Color',1.1*[0 0.447 0.741],'Linewidth',3);
set(gca,'xgrid','on','ygrid','on','tickdir','out','box','off')


%% Export movie to AVI

% AnimExport = VideoWriter('Heat1DCase3-animation.mp4','MPEG-4');
% AnimExport.Quality = 100;

% AnimExport = VideoWriter('Heat1DCase3-animation.avi','Uncompressed AVI');

AnimExport.FrameRate = 15;
open(AnimExport);
writeVideo(AnimExport,MovAnim);
close(AnimExport);
