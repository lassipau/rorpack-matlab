%% A boundary controlled 2D wave equation on an annulus
% The example case from the paper "Approximate robust output regulation of
% boundary control systems" by J.-P. Humaloja, M. Kurula and L. Paunonen.
% based on the code from https://codeocean.com/capsule/5706891/tree/v1

Nvals = 8;
Mvals = 12; % from 0 to Mvals-1

% Construction of the 2D wave equation model
[Sys,svals,sfvals,phinm1,phinm2,phiRnm,psiTHm1,psiTHm2,Q] = ConstrWave2DCase1(Nvals,Mvals);
P = @(s) Sys.C*((s*eye(2*Nvals*(2*Mvals-1))-Sys.A)\Sys.B); % transfer function approx.

%% exosystem, controller, closed-loop system
% spatial parts of the reference signal
yth1 = @(theta) -0.5/pi^2*(theta-pi).^2;
yth2 = @(theta) -0.5*sin(theta/2);
% basis for the output space
psi1 = @(m,theta) (m == 0)/sqrt(2*pi) + (m>0)/sqrt(pi).*cos(m.*theta);
psi2 = @(m,theta) 1/sqrt(pi)*sin(m.*theta);

% projection to the output space basis
dimU = 2*Mvals-1; % = dim(Y)
th = linspace(0,2*pi,100);
dth = 2*pi/numel(th);
Fc1 = zeros(dimU,1); 
Fc2 = Fc1;
for k = 0:(Mvals-1)
  if k == 0
    Fc1(1) = dth*trapz(psi1(k,th).*yth1(th));
    Fc2(1) = dth*trapz(psi1(k,th).*yth2(th));
  else
    Fc1(2*(k-1)+2:2*k+1) = ...
        dth*[trapz(psi1(k,th).*yth1(th)); trapz(psi2(k,th).*yth1(th))];
    Fc2(2*(k-1)+2:2*k+1) = ...
        dth*[trapz(psi1(k,th).*yth2(th)); trapz(psi2(k,th).*yth2(th))];
  end
end

% exosystem
ws = [-2 -1 1 2]*pi*1i; dimW = numel(ws);
S = diag(ws);
E = 0.5*sqrt(pi)*[0; 1; 0; zeros(dimU-3,1)]*[-1/2i 0 0 1/2i] + ...
    0.5*sqrt(pi)*[0; 0; 1; zeros(dimU-3,1)]*[0 1/2 1/2 0];
F = -Fc1*[0 -1/2i 1/2i 0] - Fc2*[1/2 0 0 1/2];
E = E - Q*F; % this is E_s

% controller, feedthrough Q is chosen when the plant is constructed
epsilon = 0.15; % tuning parameter
dimY = 11; % this is dim(Y_N)
dimZ = dimY*dimW;
PN = eye(dimY,dimU);
tmp = ones(dimY,1)*ws;
G1 = diag(tmp(:));
G2 = repmat(-eye(dimY)*PN, dimW, 1);
K = zeros(dimU,dimZ);
for k = 1:dimW
    K(:, (k-1)*dimY+1:k*dimY) = epsilon*(P(ws(k))\(PN'/(PN*PN')));
end

% closed-loop system
Ae = [Sys.A, Sys.B*K; G2*Sys.C, G1];
Be = [Sys.B*E; G2*F];
Ce = [Sys.C, zeros(dimU, dimZ)];
De = F;
if max(real(eig(Ae))) >= 0
    error('The closed-loop system is unstable.')
end

%% simulation
% general initial conditions for the plant
% initial conditions as functions
x0Rfun = @(r) 2*(1-(r-2).^2);
x0THfun = @(theta) (2*pi)^(-2)*theta.*(2*pi-theta).^2;
xd0Rfun = @(r) zeros(size(r));
xd0THfun = @(theta) zeros(size(theta));
x0fun = @(r,theta) x0Rfun(r).*x0THfun(theta);
xd0fun = @(r,theta) xd0Rfun(r).*xd0THfun(theta);

% projection to the state space basis
x0Rcoeffs = zeros(1,Nvals*(2*Mvals-1));
x0THcoeffs = zeros(1,Nvals*(2*Mvals-1));
xd0Rcoeffs = zeros(1,Nvals*(2*Mvals-1));
xd0THcoeffs = zeros(1,Nvals*(2*Mvals-1));
for indn = 1:Nvals
  for indm = 0:(Mvals-1)
    if indm==0
      indran = (indn-1)*(2*Mvals-1) + 1;
      x0Rcoeffs(1,indran) = ...
          integral(@(r) x0Rfun(r).*conj(phiRnm(r,indn,0)).*r,1,2);
      xd0Rcoeffs(1,indran) = ...
          integral(@(r) xd0Rfun(r).*conj(phiRnm(r,indn,0)).*r,1,2);
      
      if indn==1
        x0THcoeffs(1,indran) = integral(@(theta) ...
            x0THfun(theta).*conj(psiTHm1(theta,0)),0,2*pi);
        xd0THcoeffs(1,indran) = integral(@(theta) ...
            xd0THfun(theta).*conj(psiTHm1(theta,0)),0,2*pi);
      else
        x0THcoeffs(1,indran) = x0THcoeffs(1,1);
        xd0THcoeffs(1,indran) = xd0THcoeffs(1,1);
      end
    else
      indran = (indn-1)*(2*Mvals-1) + 2*indm - 1 + (1:2);
      x0Rcoeffs(1,indran) ...
          = [1 1]*integral(@(r) x0Rfun(r).*phiRnm(r,indn,indm).*r,1,2);
      xd0Rcoeffs(1,indran) ... 
          = [1 1]*integral(@(r) xd0Rfun(r).*phiRnm(r,indn,indm).*r,1,2);
      
      if indn==1
        x0THcoeffs(1,indran) = [integral(@(theta) x0THfun(theta).* ...
            conj(psiTHm1(theta,indm)),0,2*pi),integral(@(theta) ...
            x0THfun(theta).*conj(psiTHm2(theta,indm)),0,2*pi)];
        xd0THcoeffs(1,indran) = [integral(@(theta) xd0THfun(theta).* ...
            conj(psiTHm1(theta,indm)),0,2*pi),integral(@(theta)  ...
            xd0THfun(theta).*conj(psiTHm2(theta,indm)),0,2*pi)];
      else
        x0THcoeffs(1,indran) = x0THcoeffs(1,2*indm+(0:1));
        xd0THcoeffs(1,indran) = xd0THcoeffs(1,2*indm+(0:1));
      end
    end
  end
end
x0 = [x0Rcoeffs.*x0THcoeffs;xd0Rcoeffs.*xd0THcoeffs];
x0 = x0(:);
z0 = zeros(dimZ, 1);
v0 = ones(dimW, 1);
% to recreate the simulation in the paper, uncomment the line below
% x0 = zeros(2*Nvals*(2*Mvals-1), 1);
xe0 = [x0; z0];

v = @(t) diag(exp(ws.*t))*v0; % exosystem simulation
% solution
tspan =[0 21];
odefun = @(t,x) Ae*x + Be*v(t);
options = odeset('Jacobian', Ae);
sol = ode15s(odefun, tspan, xe0, options);

%% regulation error, output profile
th2 = linspace(0,2*pi,50);
Nt = 30;
tt = linspace(tspan(1),tspan(2),Nt*tspan(2));
xe = deval(sol,tt);
e = zeros(dimU, numel(tt));
r = e;
en = zeros(1, numel(tt));
ent = zeros(1, numel(tt)-Nt);
y = zeros(dimU, numel(tt));
d = zeros(numel(th2), numel(tt));
ref = zeros(numel(th2), numel(tt));
for k = 1:numel(tt)
  ref(:,k) = real((yth1(th2)'*[0 -1/2i 1/2i 0] + ...
      +yth2(th2)'*[1/2 0 0 1/2])*v(tt(k)));
  y(:,k) = Ce*xe(:,k);
  d(:,k) = real(cos(th2)'*[-1/2i 0 0 1/2i]*v(tt(k)) ...
      + sin(th2)'*[0 1/2 1/2 0]*v(tt(k)));
  r(:,k) = -F*v(tt(k));
end
% output profile
Y = zeros(numel(th2), numel(tt));
for s = 1:numel(tt)
  for k = 0:(Mvals-1)
    if k == 0
      Y(:,s) = Y(:,s) + y(k+1,s)*psi1(k,th2).';
    else
      Y(:,s) = Y(:,s) + y(2*k,s)*psi1(k,th2).' + ...
          y(2*k+1,s)*psi2(k,th2).';
    end
  end
  en(s) = 2*pi/numel(th2)*norm(real(Y(:,s)-ref(:,s)))^2;
end
% time average of the regulation error
for s = 1:numel(tt)-Nt
  ent(s) = trapz(0:1/Nt:1, en(s:s+Nt));
end
Y = real(Y); % neglect imaginary part (should be zero)

%% variables for visualization
% annular grid
Nr = 50;
Nth = 50;
[rrg,thg] = meshgrid(linspace(1,2,Nr),linspace(0,2*pi,Nth));
rr = rrg(:).';
th = thg(:).';
xxg = reshape(rr.*cos(th),Nth,Nr);
yyg = reshape(rr.*sin(th),Nth,Nr);

% basis for the wave profile
pnvals = zeros(Nvals*(2*Mvals-1),numel(rr));
for indn = 1:Nvals
  for indm = 0:(Mvals-1)
    if indm==0
      indran = (indn-1)*(2*Mvals-1) + 1;
      pnvals(indran,:) = phinm1(rr,th,indn,0);
    else
      indran = (indn-1)*(2*Mvals-1) + 2*indm - 1 + (1:2);
      pnvals(indran,:) = [phinm1(rr,th,indn,indm);phinm2(rr,th,indn,indm)];
    end
  end
end

%% visualizations
% regulation error
figure(1)
del = 0.01; % approximate tracking tolerance
semilogy(tt(1:numel(tt)-Nt), real(ent), 'linewidth', 2)
xlim([0, tspan(2)-1])
hold on
semilogy([0, tspan(2)-1], del*norm(v0)^2*[1 1],'--k', 'linewidth', 1.5)
xlabel('$t$', 'interpreter', 'latex', 'fontsize', 13)
ylabel('Tracking errors $\int_t^{t+1}\|e(s)\|^2ds$', ...
    'interpreter', 'latex', 'fontsize', 13)
hold off
title('The regulation error.','interpreter', 'latex', 'fontsize', 15)
% print -dpng /results/error

% output profile
figure(2)
surf(tt, th2, Y)
xlabel('$t$', 'interpreter', 'latex', 'fontsize', 13)
ylabel('$\theta$', 'interpreter', 'latex', 'fontsize', 13)
zlabel('$y$', 'interpreter', 'latex', 'fontsize', 13, 'rotation', 0)
xlim([0, 10])
title('The output profile.','interpreter', 'latex', 'fontsize', 15)
% print -dpng /results/output

% reference profile
figure(3)
surf(tt, th2, ref)
xlabel('$t$', 'interpreter', 'latex', 'fontsize', 13)
ylabel('$\theta$', 'interpreter', 'latex', 'fontsize', 13)
zlabel('$y_{ref}$', 'interpreter', 'latex', 'fontsize', 13, 'rotation', 0)
xlim([0, 10])
title('The reference profile.','interpreter', 'latex', 'fontsize', 15)
% print -dpng /results/reference

% distubance profile
figure(4)
surf(tt, th2, d)
xlabel('$t$', 'interpreter', 'latex', 'fontsize', 13)
ylabel('$\theta$', 'interpreter', 'latex', 'fontsize', 13)
zlabel('$d$', 'interpreter', 'latex', 'fontsize', 13, 'rotation', 0)
xlim([0 6])
title('The disturbance profile.','interpreter', 'latex', 'fontsize', 15)
% print -dpng /results/disturbance

% state at a specific time t
figure(5)
ind = 9*Nt+1;
ind = 10*Nt+1;
tt(ind)
zzvalmat = xe(1:2:2*Nvals*(2*Mvals-1),ind).'*pnvals;
surf(xxg,yyg,reshape(real(zzvalmat),Nth,Nr));
xlabel('$\zeta_1$', 'interpreter', 'latex', 'fontsize', 13)
ylabel('$\zeta_2$', 'interpreter', 'latex', 'fontsize', 13)
zlabel('$w$', 'interpreter', 'latex', 'fontsize', 13, 'rotation', 0)
title(['The state of the system at $t=' num2str(round(tt(ind),1)) '$.'],'interpreter', 'latex', 'fontsize', 15);
% print -dpng /results/state9

%% Animation
% compute the wave profile
zzvalmat = zeros(numel(tt),Nth*Nr);
for ind = 1:numel(tt)
  zzvalmat(ind,:) = xe(1:2:2*Nvals*(2*Mvals-1),ind).'*pnvals;
end
zzvalmat = real(zzvalmat);

% Animate
figure(6)
zlims =  [min(zzvalmat(:)) max(zzvalmat(:))];
axlims = [-2 2 -2 2 zlims];
obj = VideoWriter('wave.avi');
set(obj, 'FrameRate', 25, 'Quality', 25);
open(obj)
for ind = 1:2:numel(tt)
  clf
  surf(xxg,yyg,reshape(zzvalmat(ind,:),Nth,Nr));
  title(['t=', num2str(tt(ind), ...
      2+double(tt(ind)>=1)+max(0,floor(log10(tt(ind)))))])
  axis(axlims)
  caxis(zlims)
  xlabel('$x$','Interpreter','latex','Fontsize',13)
  ylabel('$y$','Interpreter','latex','Fontsize',13)
  frame = getframe(gcf);
  writeVideo(obj, frame)
end
close(obj)