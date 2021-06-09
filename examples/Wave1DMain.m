%% A Wave equation in 1D with Neumann boundary control and Dirichlet observation, 
% based on Guo & Krstic IFAC 2017, Main file
% Simulation based on modal approximation

% Unstable system, noncollocated I/O. Can be stabilized with output
% feedback, but the low-gain controller achieves a poor stability margin

addpath(genpath('../RORPack/'))

N = 60; 

% Initial state of the plant
w0fun = @(x) zeros(size(x));
%w0fun = @(x) 1+cos(3*pi*x)+cos(6*x);
wd0fun = @(x) zeros(size(x));
%x0fun = @(x,y) 0.5*(1+cos(pi*(1-x))).*(1-1/4*cos(2*pi*y));
% x0fun = @(x,y) 0.5*(1+cos(pi*(1-x)));
%x0fun = @(x,y) 1/2*x.^2.*(3-2*x)-1;
%x0fun = @(x,y) 1/2*x.^2.*(3-2*x)-1/2;
%x0fun = @(x,y) 1*(1-x).^2.*(3-2*(1-x))-1;
%x0fun = @(x,y) .5*(1-x).^2.*(3-2*(1-x))-.5;
%x0fun = @(x,y) 1/4*(x.^3-1.5*x.^2)-1/4;
%x0fun = @(x,y) .2*x.^2.*(3-2*x)-.5;

[~,Sys,phin,Kinf,Linf] = Constr1DWave(w0fun,wd0fun,N);

%yref = @(t) sin(2*t)+.1*cos(6*t);
%yref = @(t) sin(2*t)+.2*cos(3*t);
%yref = @(t) ones(size(t));

%wdist = @(t) ones(size(t));
%wdist = @(t) sin(t);
wdist = @(t) zeros(size(t));

% Triangle signal case
% Begin by defining the function on a single period 0<t<2
% A nonsmooth triangle signal
yref1per = @(t) (2.*t-1).*(t>=0).*(t<=1)+(3-2.*t).*(t>1).*(t<2);
% The constant part of the signal cannot be tracked due to the second order
% zero of the plant at zero.
% We therefore normalize yref(t) to have average zero on [0,2]
yr_ave = integral(yref1per,0,2);
yref = @(t) yref1per(mod(t,2)) - yr_ave/2;

% Case 1:
% yref = @(t) .5*sin(2*t)+.5*cos(3*t);
%yref = @(t) sin(pi/2*t)-2*cos(pi/2*t);
% yref = @(t) zeros(size(t));
%wdist = @(t) zeros(size(t));
%wdist = @(t) sin(2*pi*t);

% Case 2:
% yref = @(t) ones(size(t));
% wdist = @(t) ones(size(t));

% Case 3:
% yref = @(t) sin(2*t)+.1*cos(6*t);
% wdist = @(t) sin(t);


% freqs = [-3i -2i -1i 0 1i 2i 3i];

% freqs = [1 2 3 4];
freqs = [pi/2,pi,2*pi,3*pi,4*pi];

% dimX = size(Sys.A,1);
% Pappr = @(s) Sys.C*((s*eye(dimX)-Sys.A)\Sys.B)+Sys.D;
% 
% Pvals = cell(1,length(freqs));
% for ind = 1:length(freqs)
%   Pvals{ind} = Pappr(freqs(ind));
% end

% for ind = 1:length(freqsReal)
%   if cond(Pappr(freqs(ind)))>1e6
%     warning(['The matrix P(iw_k) for k=' num2str(ind) ' is nearly singular!'])
%   end
% end


% epsgainrange = [10,50];
% epsgain = 13;
% [ContrSys,epsgain] = ConstrContrLG(freqs,Pvals,epsgain,Sys);
%[ContrSys,epsgain] = ConstrContrPassive(freqs,Pvals,epsgainrange,Sys);

% % Step 1: Compute stabilizing operators K2 and L so that A+B*K2 and A+L*C
% % are exponentially stable
% K21 = -lqr(full(Sys.A),Sys.B,eye(size(Sys.A)),1*eye(dimU),zeros(dimX,dimU));
% L = -lqr(full((Sys.A).'),(Sys.C).',eye(size(Sys.A)),0.1*eye(dimY),zeros(dimX,dimY)).';

% K0 = -1/3*[3,2,zeros(1,2*N-2)];
% L0 = -.6*[3;2;zeros(2*N-2,1)];
% L0 = -.3*[0;1;zeros(2*N-2,1)];

% Prestabilization
k_m = .3;
Sys.A = Sys.A-k_m*Sys.B*Sys.Cm;

% CHEATING: Prestabilize the Jordan block at s=0
% Sys.A(1:2,1:2) = [-1 1;0 0];

kappa = .9;
ell = .8;
% K_S = -kappa*(Sys.B)'+K0;
% K_S = -kappa*Kinf+K0;
K_S = -kappa*Kinf;
% PlotEigs(Sys.A+Sys.B*K_S,[-1.5, 0.1, NaN, NaN])
 

% L = [zeros(dimX,1),-ell*(Sys.Cm)'+L0];
% L = [L0,-ell*(Sys.Cm)'];
L = -ell*Linf;
% L = L0;

% PlotEigs(Sys.A+L*[Sys.C;Sys.Cm],[-1 0 -10 10])
% PlotEigs(Sys.A+L*[Sys.C;Sys.Cm])
% PlotEigs(Sys.A+L*Sys.C,[-1.5, .01, NaN, NaN])
% PlotEigs(Sys.A+L*Sys.C)



% [ContrSys,K21] = ConstrContrObsBased(freqs,Sys,K_S,L,'poleplacement',3);
[ContrSys,K21] = ConstrContrObsBased(freqs,Sys,K_S,L,'LQR',3);
% [ContrSys,K21] = ConstrContrDualObsBased(freqs,Sys,K_S,L,'LQR',3);
% [ContrSys,K21] = ConstrContrDualObsBased(freqs,Sys,K_S,L,'poleplacement',3);

CLSys = ConstrCLSys(Sys,ContrSys);

% PlotEigs(full(CLSys.Ae),[-2 NaN NaN NaN])
stabmarg = CLStabMargin(CLSys)

%%

% w0fun = @(x) zeros(size(x));
% w0fun = @(x) 1+cos(3*pi*x)+cos(6*x);
% w0fun = @(x) -cos(2*pi*x);
w0fun = @(x) 30*x.^2.*(1-x).^2-1;
wd0fun = @(x) zeros(size(x));
%x0fun = @(x,y) 0.5*(1+cos(pi*(1-x))).*(1-1/4*cos(2*pi*y));
% x0fun = @(x,y) 0.5*(1+cos(pi*(1-x)));
%x0fun = @(x,y) 1/2*x.^2.*(3-2*x)-1;
%x0fun = @(x,y) 1/2*x.^2.*(3-2*x)-1/2;
%x0fun = @(x,y) 1*(1-x).^2.*(3-2*(1-x))-1;
%x0fun = @(x,y) .5*(1-x).^2.*(3-2*(1-x))-.5;
%x0fun = @(x,y) 1/4*(x.^3-1.5*x.^2)-1/4;
%x0fun = @(x,y) .2*x.^2.*(3-2*x)-.5;

x0 = Constr1DWave(w0fun,wd0fun,N);

xe0 = [x0;zeros(size(ContrSys.G1,1),1)];

Tend = 20;
tgrid = linspace(0,Tend,200);



CLsim = SimCLSys(CLSys,xe0,yref,wdist,tgrid,[]);

% Choose whther or not to print titles of the figures
PrintFigureTitles = true;

figure(1)
subplot(2,1,1)
plotOutput(tgrid,yref,CLsim,PrintFigureTitles)
subplot(2,1,2)
plotErrorNorm(tgrid,CLsim,PrintFigureTitles)

q = length(freqs);
if freqs(1)==0,dimZ = 2*q-1; else dimZ=2*q;end
dimX = size(Sys.A,1);
dimU = size(ContrSys.K,1);

K1full = [ContrSys.K(:,1:dimZ), zeros(dimU,dimX)];
K2full = [zeros(dimU,dimZ) ContrSys.K(:,(dimZ+1):end)];
figure(2)
plot(tgrid,[ContrSys.K*CLsim.xesol((2*N+1):end,:);K1full*CLsim.xesol((2*N+1):end,:);K2full*CLsim.xesol((2*N+1):end,:)],'Linewidth',2);
obserror = sum((CLsim.xesol(1:2:(2*N),:)-CLsim.xesol(dimX+dimZ+(1:2:(2*N)),:)).^2,1);
figure(3)
plot(tgrid,obserror,'Linewidth',2);

figure(4)
colormap jet
%PlotHeat2DSurf(x0,spgrid,[-1.4,1.4])
spgrid = linspace(0,1,N);
Plot1DWaveSurf(CLsim.xesol(1:2*N,:),phin,spgrid,tgrid)
% Plot1DWaveSurf(CLsim.xesol(1:2*N,:)-CLsim.xesol(dimX+dimZ+(1:(2*N)),:),phin,spgrid,tgrid,[-9 9])
% Plot1DWaveSurf(CLsim.xesol(dimX+dimZ+(1:(2*N)),:),phin,spgrid,tgrid,[-4 4])
% set(gca,'ztick',-8:4:8);
%colormap jet

% figure(5)
% tt = linspace(0,16,500)
% plot(tt,yref(tt),'Color',1.1*[0 0.447 0.741],'Linewidth',3);
% set(gca,'xgrid','on','ygrid','on','tickdir','out','box','off')

%% Animation


% figure(4)
% colormap jet
% No movie recording
% [~,zlims] = Anim1DWaveSpectral(CLsim.xesol(1:2*N,:),phin,spgrid,tgrid,0.03,0);

% Movie recording
% [MovAnim,zlims] = AnimHeat2Dtest1(CLsim,spgrid,tgrid,0,1);

% Tpause = 0.05;
% record = 0;
% [MovAnim,zlims] = Anim1DWaveSpectral(CLsim.xesol(1:2*N,:),phin,spgrid,tgrid,Tpause,record)

%movie(MovAnim)

% figure(3)
% colormap jet
% %PlotHeat2DSurf(x0,spgrid,[-1.4,1.4])
% PlotHeat2DSurf(x0,spgrid,zlims)
% 
% figure(4)
% tt = linspace(0,16,500)
% plot(tt,yref(tt),'Color',1.1*[0 0.447 0.741],'Linewidth',3);
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


%%

% K21=K21(:);
% phinvals = phin(spgrid,(0:(N-1)).'); 
% K21(1:2:end)*phinvals;
% plot(spgrid,K21(2:2:end)*phinvals);


% A = Sys.A;
% B = Sys.B;
% C = Sys.C;
% C2 = Sys.C2;

% max(max(real(eig(full(A+B*K21)))))
% max(max(real(eig(full(A-0.75*B*(C2+C))))))
% max(max(real(eig(full(A-75*B*C)))))

