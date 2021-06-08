%% A Wave equation in 1D with Neumann boundary control and Dirichlet observation, 
% based on Guo & Krstic IFAC 2017, Main file
% Simulation based on modal approximation

% Unstable system, noncollocated I/O. 


N = 200; 

% Initial state of the plant
w0fun = @(x) zeros(size(x));
% w0fun = @(x) 1+cos(3*pi*x)-cos(6*x);
wd0fun = @(x) zeros(size(x));
%x0fun = @(x,y) 0.5*(1+cos(pi*(1-x))).*(1-1/4*cos(2*pi*y));
% x0fun = @(x,y) 0.5*(1+cos(pi*(1-x)));
%x0fun = @(x,y) 1/2*x.^2.*(3-2*x)-1;
%x0fun = @(x,y) 1/2*x.^2.*(3-2*x)-1/2;
%x0fun = @(x,y) 1*(1-x).^2.*(3-2*(1-x))-1;
%x0fun = @(x,y) .5*(1-x).^2.*(3-2*(1-x))-.5;
%x0fun = @(x,y) 1/4*(x.^3-1.5*x.^2)-1/4;
%x0fun = @(x,y) .2*x.^2.*(3-2*x)-.5;
tic
% The example uses the same constructor as another Wave example, but now
% with only one output (we only use the velocity measurement, for simplicity)
[x0,Sys,phin] = Constr1DWave2(w0fun,wd0fun,N);
A = Sys.A;
toc
% % Debugging
% B = Sys.B;
% C = Sys.C;
% Cm = Sys.Cm;
% D = Sys.D;
% 
% Am = A-k_m*B*Cm;
% 
% Resolv = @(s,T) inv(s*eye(size(T))-T);
% 
% PKappr = @(s) (C+D*K_S)*((s*speye(2*N)-(Am+B*K_S))\B)+D;
% gammaappr = @(s) PKappr(s)^(-1)*((C+D*K_S)*((s*eye(2*N)-Am-B*K_S)\L)-1);
% 
% CRB = @(s) C*((s*eye(2*N)-A)\B);
% CRB = @(s) 2/(exp(s)-exp(-s));
% 
% CmRB = @(s) Cm*((s*eye(2*N)-A)\B);
% CmRB = @(s) 1/s^2;
% 
% KSRB = @(s) K_S*((s*eye(2*N)-A)\B);
% KSRB = @(s) -kappa*(1+exp(-2*s))/(1-exp(-2*s));
% 
% Pm = @(s) C*((s*eye(2*N)-Am)\B);
% Pm = @(s) C*((s*eye(2*N)-A)\B)-k_m*C*((s*eye(2*N)-A)\B)*(1+k_m*Cm*((s*eye(2*N)-A)\B))^(-1)*Cm*((s*eye(2*N)-A)\B);
% Pm = @(s) CRB(s)-k_m*CRB(s)*(1+k_m*CmRB(s))^(-1)*CmRB(s);
% 
% 
% KSRmB = @(s) K_S*((s*eye(2*N)-Am)\B);
% KSRmB = @(s) K_S*((s*eye(2*N)-A)\B)-k_m*K_S*((s*eye(2*N)-A)\B)*(1+k_m*Cm*((s*eye(2*N)-A)\B))^(-1)*Cm*((s*eye(2*N)-A)\B);
% KSRmB = @(s) KSRB(s)-k_m*KSRB(s)*(1+k_m*CmRB(s))^(-1)*CmRB(s);
% 
% CRmL = @(s) C*((s*eye(2*N)-Am)\L);
% CRmL = @(s) C*((s*eye(2*N)-A)\L)-k_m*C*((s*eye(2*N)-A)\B)*(1+k_m*Cm*((s*eye(2*N)-A)\B))^(-1)*Cm*((s*eye(2*N)-A)\L);
% 
% 
% KSRmL = @(s) K_S*((s*eye(2*N)-Am)\L);
% KSRmL = @(s) K_S*((s*eye(2*N)-A)\L)-k_m*K_S*((s*eye(2*N)-A)\B)*(1+k_m*Cm*((s*eye(2*N)-A)\B))^(-1)*Cm*((s*eye(2*N)-A)\L);
% 
% gammaappr2 = @(s) KSRmL(s)+(KSRmB(s)-1)*Pm(s)^(-1)*(1-CRmL(s));
% 
% % PlotEigs(A-.3*(B*B'),[-3 0.1 -30 30])
% % Cm*((s*eye(2*N)-A)\B)
% % Cm*((s*eye(2*N)-A)\L)
% % C*((s*eye(2*N)-A)\B)
% % C*((s*eye(2*N)-A)\L)
% 
% ss = linspace(1,2);
% vals = zeros(size(ss));
% for ind = 1:length(ss)
% %   vals(ind)=abs(C*((1i*ss(ind)*eye(2*N)-A)\B)-2/(exp(1i*ss(ind))-exp(-1i*ss(ind))));
% % vals(ind)=abs(K_S*((1i*ss(ind)*eye(2*N)-A)\B)-(-kappa)*(1+exp(-2*1i*ss(ind)))/(1-exp(-2*1i*ss(ind))));
% %   vals(ind)=abs(C*((1i*ss(ind)*eye(2*N)-A)\L)-(-ell)*(1+exp(-2*1i*ss(ind)))/(1-exp(-2*1i*ss(ind))));
% %   vals(ind)=abs(K_S*((1i*ss(ind)*eye(2*N)-A)\L)-(kappa*ell)*2/(exp(1i*ss(ind))-exp(-1i*ss(ind))));
% 
% %    vals(ind)=abs(Cm*((1i*ss(ind)*eye(2*N)-A)\B)-1/((1i*ss(ind))^2));
% 
% %   vals(ind)=abs(C*((1i*ss(ind)*eye(2*N)-Am)\B)-1i*ss(ind)^2/(sin(ss(ind))*(-ss(ind)^2+k_m)));
%   
% %   vals(ind)=abs(K_S*((1i*ss(ind)*eye(2*N)-Am)\B)- kappa*(1+exp(-(2*1i)*ss(ind)))*ss(ind)^2/((-1+exp(-(2*1i)*ss(ind)))*(ss(ind)^2-k_m)));
% 
% %   vals(ind)=abs(C*((1i*ss(ind)*eye(2*N)-Am)\L)-(1/2*1i)*((ss(ind)^2-k_m)*exp(-(3*1i)*ss(ind))+2*k_m*exp(-(2*1i)*ss(ind))+(-ss(ind)^2+k_m)*exp(1i*ss(ind))-2*k_m)*ell/((-1+exp(-(2*1i)*ss(ind)))*sin(ss(ind))*(ss(ind)^2-k_m)));
% 
% %   vals(ind)=abs(K_S*((1i*ss(ind)*eye(2*N)-Am)\L)-(-1i)*kappa*ell*((ss(ind)^2-k_m)*exp(-(2*1i)*ss(ind))-ss(ind)^2+(1/2)*exp(-(3*1i)*ss(ind))*k_m-(1/2)*exp(1i*ss(ind))*k_m+k_m)/(sin(ss(ind))*(-1+exp(-(2*1i)*ss(ind)))*(ss(ind)^2-k_m)));
% 
% % vals(ind) = norm(Resolv(1i*ss(ind),Am)-(Resolv(1i*ss(ind),A)-k_m*Resolv(1i*ss(ind),A)*B*(1+k_m*Cm*Resolv(1i*ss(ind),A)*B)^(-1)*Cm*Resolv(1i*ss(ind),A)));
%  vals(ind) = abs(gammaappr(1i*ss(ind))-gammaappr2(1i*ss(ind)));
% end
% plot(ss,vals)

%




%yref = @(t) sin(2*t)+.1*cos(6*t);
%yref = @(t) sin(2*t)+.2*cos(3*t);
%yref = @(t) ones(size(t));

%wdist = @(t) ones(size(t));
%wdist = @(t) sin(t);
%wdist = @(t) zeros(size(t));

% Case 1:
%yref = @(t) [sin(2*t)+.5*cos(1*t);zeros(size(t))]
yref = @(t) 3*sin(1*t)+4*cos(2*t);
% yref = @(t) sin(pi/4*t)-cos(pi/4*t);
% yref = @(t) zeros(size(t));
wdist = @(t) zeros(size(t));
% wdist = @(t) 3*sin(4*t);
% wdist = @(t) .5*cos(7*t)+.5*cos(8*t);

% Case 2:
% yref = @(t) ones(size(t));
% wdist = @(t) ones(size(t));

% Case 3:
% yref = @(t) sin(2*t)+.1*cos(6*t);
% wdist = @(t) sin(t);


% freqs = [-3i -2i -1i 0 1i 2i 3i];
% freqsReal = [0 1 2 3 4];
freqsReal = [1 2 7 8];
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

% We choose the stabilizing operators as K_S = -kappa*(B'+K0), L = -ell*(C+K0.')' with
% kappa,ell>0 and K0 only affects the eigenspace corresponding to w_0=0
% In this example, gamma(s) can be computed explicitly as
% gamma(iw) = -(kappa+ell)*cos(w) + 1i*(-sin(w)*(kappa*ell+1))


% K0 = -1/3*[3,2,zeros(1,2*N-2)];
% L0 = -.6*[3;2;zeros(2*N-2,1)];
% L0 = -.3*[0;1;zeros(2*N-2,1)];

k_m = .25;
% Sys.A = A-k_m*Sys.B*Sys.Cm;


kappa = .9;
ell = .8;
% K_S = -kappa*(Sys.B)'+K0;
K_S = -kappa*(Sys.B)';
% PlotEigs(Sys.A+Sys.B*K_S)


% L = [zeros(dimX,1),-ell*(Sys.Cm)'+L0];
% L = [L0,-ell*(Sys.Cm)'];
L = -ell*(Sys.C)';

% PlotEigs(Sys.A+L*[Sys.C;Sys.Cm],[-1 0 -10 10])
%  PlotEigs(Sys.A+L*[Sys.C;Sys.Cm])
PlotEigs(Sys.A+L*Sys.C,[-.4, .01, NaN, NaN])
% PlotEigs(Sys.A+L*Sys.C)

%%

%gammafun = @(w) -(kappa+ell)*cos(w) + 1i*(-sin(w)*(kappa*ell+1));
gammafun = @(w) (cos(w).*(-(ell+kappa).*w.^2+k_m*ell)-k_m*ell)./w.^2 + 1i*sin(w).*(k_m-(ell*kappa+1).*w.^2)./w.^2;

tic
[ContrSys,K_S] = ConstrContrREObsReal(freqsReal,Sys,K_S,L,gammafun);
% [ContrSys,K_S] = ConstrContrREObsReal(freqsReal,Sys,K_S,L); % gammafun obtained using the numerical approximation
toc
%
% eig(full(ContrSys.G1))

% epsgain

% Sys.C = [Sys.C;Sys.Cm];
% Sys.D = [Sys.D;Sys.Dm];
tic
CLSys = ConstrCLSys(Sys,ContrSys);
toc
stabmarg = CLStabMargin(CLSys)

%%
xe0 = [x0;zeros(size(ContrSys.G1,1),1)];

Tend = 28;
tgrid = linspace(0,Tend,300);


tic
CLsim = SimCLSys(CLSys,xe0,yref,wdist,tgrid,[]);
toc


figure(1)
plotBasics(tgrid,yref,CLsim)

% figure(2)
% q = length(freqs);
% if freqs(1)==0,dimZ = 2*q-1; else dimZ=2*q;end
% dimX = size(Sys.A,1);
% dimU = size(ContrSys.K,1);
% 
% K1full = [ContrSys.K(:,1:dimZ), zeros(dimU,dimX)];
% K2full = [zeros(dimU,dimZ) ContrSys.K(:,(dimZ+1):end)];
% plot(tgrid,[ContrSys.K*CLsim.xesol((2*N+1):end,:);K1full*CLsim.xesol((2*N+1):end,:);K2full*CLsim.xesol((2*N+1):end,:)],'Linewidth',2);
% obserror = sum((CLsim.xesol(1:2:(2*N),:)-CLsim.xesol(dimX+dimZ+(1:2:(2*N)),:)).^2,1);
% plot(tgrid,obserror,'Linewidth',2);

figure(3)
colormap jet
%PlotHeat2DSurf(x0,spgrid,[-1.4,1.4])
spgrid = linspace(0,1,100);

%  Plot1DWaveSurf(CLsim.xesol(1:2*N,:),phin,spgrid,tgrid)
% Less points in the time-axis
Plot1DWaveSurf(CLsim.xesol(1:2*N,1:2:end),phin,spgrid(1:2:end),tgrid(1:2:end))
% Plot1DWaveSurf(CLsim.xesol(1:2*N,:)-CLsim.xesol(dimX+dimZ+(1:(2*N)),:),phin,spgrid,tgrid,[-9 9])
% Plot1DWaveSurf(CLsim.xesol(dimX+dimZ+(1:(2*N)),:),phin,spgrid,tgrid,[-4 4])
set(gca,'ztick',-8:4:8);
%colormap jet

% figure(4)
% tt = linspace(0,16,500)
% plot(tt,yref(tt),'Color',1.1*[0 0.447 0.741],'Linewidth',3);
% set(gca,'xgrid','on','ygrid','on','tickdir','out','box','off')

%%


% figure(2)
% colormap jet
% No movie recording
% [~,zlims] = AnimHeat2Dtest1(CLsim,spgrid,tgrid,0.03,0);

% Movie recording
% [MovAnim,zlims] = AnimHeat2Dtest1(CLsim,spgrid,tgrid,0,1);

%movie(MovAnim)

figure(3)
colormap jet
PlotHeat2DSurf(x0,spgrid,[-1.4,1.4])
%PlotHeat2DSurf(x0,spgrid,zlims)

figure(4)
tt = linspace(0,16,500);
plot(tt,yref(tt),'Color',1.1*[0 0.447 0.741],'Linewidth',3);
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

