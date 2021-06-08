%% ACC 2019 paper simulations with Sepideh

addpath(genpath('../RORPack/'))

load('ACCmat.mat')




N = size(A,1);



Sys.A = A;
Sys.B = B;
Sys.Bd = zeros(N,1);
Sys.C = C;
Sys.D = zeros(4);

x0 = 100*ones(N,1);

% Case 1:
yref = @(t) [1;2;-1;0.1]*(sin(20*t)+.5*cos(60*t));

yr1=@(t) .005*sin(20*t)+.005*sin(60*t);
yr2=@(t) .005*cos(20*t)+.005*cos(60*t);
yr3=@(t) 7.1e-4*(1+sign(yr1(t)))/4;
yr4=@(t) 7.1e-4-yr3(t);

yref=@(t) [yr1(t);yr2(t);yr3(t);yr4(t)];

% yref = @(t) sin(pi/4*t)-cos(pi/4*t);
% yref = @(t) zeros(size(t));
wdist = @(t) zeros(size(t));

% Case 2:
% yref = @(t) ones(size(t));
% wdist = @(t) ones(size(t));

% Case 3:
% yref = @(t) sin(2*t)+.1*cos(6*t);
% wdist = @(t) sin(t);


% freqs = [-3i -2i -1i 0 1i 2i 3i];
freqsReal = [0 20 60];
% freqsReal = pi/4;
freqs = 1i*freqsReal;

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
%[ContrSys,epsgain] = ConstrContrPassiveReal(freqsReal,Pvals,epsgainrange,Sys);

% % Step 1: Compute stabilizing operators K2 and L so that A+B*K2 and A+L*C
% % are exponentially stable
% K21 = -lqr(full(Sys.A),Sys.B,eye(size(Sys.A)),1*eye(dimU),zeros(dimX,dimU));
% L = -lqr(full((Sys.A).'),(Sys.C).',eye(size(Sys.A)),0.1*eye(dimY),zeros(dimX,dimY)).';

% K0 = -1/3*[3,2,zeros(1,2*N-2)];
% L0 = -.6*[3;2;zeros(2*N-2,1)];
% % L0 = -.3*[0;1;zeros(2*N-2,1)];
% 
% % k_m = .25;
% % Sys.A = Sys.A-k_m*Sys.B*Sys.Cm;
% 
% % CHEATING: Prestabilize the Jordan block at s=0
% Sys.A(1:2,1:2) = [-1 1;0 0];
% 
% kappa = .9;
% ell = .9;
% % K_S = -kappa*(Sys.B)'+K0;
% K_S = -kappa*Kinf+K0;
% %  PlotEigs(Sys.A+Sys.B*K_S,[-1.5, 0.1, NaN, NaN])
%  
% 
% 
% 
% % L = [zeros(dimX,1),-ell*(Sys.Cm)'+L0];
% % L = [L0,-ell*(Sys.Cm)'];
% L = -ell*Linf;
% % L = L0;
% 
% % PlotEigs(Sys.A+L*[Sys.C;Sys.Cm],[-1 0 -10 10])
% %  PlotEigs(Sys.A+L*[Sys.C;Sys.Cm])
% % PlotEigs(Sys.A+L*Sys.C,[-.4, .01, NaN, NaN])
% %  PlotEigs(Sys.A+L*Sys.C)

K_S= K2;
L = L;

% [ContrSys,K21] = ConstrContrObsBasedReal(freqsReal,Sys,K_S,L,'poleplacement',3);
[ContrSys,K21] = ConstrContrObsBasedReal(freqsReal,Sys,K_S,L,'LQR',3);

% epsgain

%%
CLSys = ConstrCLSys(Sys,ContrSys);

% PlotEigs(full(CLSys.Ae),[-2 NaN NaN NaN])
stabmarg = CLStabMargin(CLSys)


xe0 = [x0;zeros(size(ContrSys.G1,1),1)];

Tend = 2;
tgrid = linspace(0,Tend,200);



CLsim = SimCLSys(CLSys,xe0,yref,wdist,tgrid,[]);

figure(1)
subplot(2,1,1)
hold off
cla
hold on
plot(tgrid,yref(tgrid),'Color',1.1*[0 0.447 0.741],'Linewidth',2);
plot(tgrid,CLsim.output,'Color', [0.85 0.325 0.098],'Linewidth',2);
title('Output $y(t)$ (red) and the reference $y_{ref}(t)$ (blue)','Interpreter','latex','Fontsize',16)
set(gca,'xgrid','off','tickdir','out','box','off')
subplot(2,1,2)
plot(tgrid,CLsim.error,'Linewidth',2);
set(gca,'xgrid','on','ygrid','on','tickdir','out','box','off')
title('Tracking  error $y(t)-y_{ref}(t)$','Interpreter','latex','Fontsize',16)
%set(gcf,'color',1/255*[252 247 255])

figure(2)
q = length(freqs);
if freqs(1)==0,dimZ = 2*q-1; else dimZ=2*q;end
dimX = size(Sys.A,1);
dimU = size(ContrSys.K,1);

K1full = [ContrSys.K(:,1:dimZ), zeros(dimU,dimX)];
K2full = [zeros(dimU,dimZ) ContrSys.K(:,(dimZ+1):end)];
plot(tgrid,ContrSys.K*CLsim.xesol((N+1):end,:),'Linewidth',2);
% obserror = sum((CLsim.xesol(1:2:(2*N),:)-CLsim.xesol(dimX+dimZ+(1:2:(2*N)),:)).^2,1);
% plot(tgrid,obserror,'Linewidth',2);

% figure(3)
% colormap jet
% %PlotHeat2DSurf(x0,spgrid,[-1.4,1.4])
% spgrid = linspace(0,1,100);
% 
% Plot1DWaveSurf(CLsim.xesol(1:2*N,:),phin,spgrid,tgrid)
% % Plot1DWaveSurf(CLsim.xesol(1:2*N,:)-CLsim.xesol(dimX+dimZ+(1:(2*N)),:),phin,spgrid,tgrid,[-9 9])
% % Plot1DWaveSurf(CLsim.xesol(dimX+dimZ+(1:(2*N)),:),phin,spgrid,tgrid,[-4 4])
% set(gca,'ztick',-8:4:8);
%colormap jet

% figure(4)
% tt = linspace(0,16,500)
% plot(tt,yref(tt),'Color',1.1*[0 0.447 0.741],'Linewidth',3);
% set(gca,'xgrid','on','ygrid','on','tickdir','out','box','off')



save('ACCcontroller.mat','Sys','ContrSys','CLSys');

%% Animation


figure(2)
colormap jet
% No movie recording
% [~,zlims] = AnimHeat2Dtest1(CLsim,spgrid,tgrid,0.03,0);

% Movie recording
% [MovAnim,zlims] = AnimHeat2Dtest1(CLsim,spgrid,tgrid,0,1);

Tpause = 0.05;
record = 0;
[MovAnim,zlims] = Anim1DWaveSpectral(CLsim.xesol(1:2*N,:),phin,spgrid,tgrid,Tpause,record)

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
phinvals = phin(spgrid,(0:(N-1)).'); 
% K21(1:2:end)*phinvals;
plot(spgrid,K21(2:2:end)*phinvals);


A = Sys.A;
B = Sys.B;
C = Sys.C;
C2 = Sys.C2;

max(max(real(eig(full(A+B*K21)))))
max(max(real(eig(full(A-0.75*B*(C2+C))))))
max(max(real(eig(full(A-75*B*C)))))

