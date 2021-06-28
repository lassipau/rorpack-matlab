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
x10fun = @(x) zeros(size(x));
x20fun = @(x) zeros(size(x));
x30fun = @(x) zeros(size(x));
x40fun = @(x) zeros(size(x));
%x0fun = @(x,y) 0.5*(1+cos(pi*(1-x))).*(1-1/4*cos(2*pi*y));
% x0fun = @(x,y) 0.5*(1+cos(pi*(1-x)));
%x0fun = @(x,y) 1/2*x.^2.*(3-2*x)-1;
%x0fun = @(x,y) 1/2*x.^2.*(3-2*x)-1/2;
%x0fun = @(x,y) 1*(1-x).^2.*(3-2*(1-x))-1;
%x0fun = @(x,y) .5*(1-x).^2.*(3-2*(1-x))-.5;
%x0fun = @(x,y) 1/4*(x.^3-1.5*x.^2)-1/4;
%x0fun = @(x,y) .2*x.^2.*(3-2*x)-.5;


[x0,spgrid,Sys] = ConstrLHMNC18(x10fun,x20fun,x30fun,x40fun,N);

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


% freqs = [-3i -2i -1i 0 1i 2i 3i];
freqsReal = [1 2];


% Check that the system, the reference and disturbance signals, and the
% frequencies are defined in a consistent way.
Sys = SysConsistent(Sys,yref,wdist,freqsReal);



%% Controller construction


% Simple passive robust controller, used in the original Timoshenko beam
% example in the LHMNC 2018 conference paper (simulation not included in
% the paper).

dimX = size(Sys.A,1);
Pappr = @(s) Sys.C*((s*eye(dimX)-Sys.A)\Sys.B)+Sys.D;

Pvals = cell(1,length(freqsReal));
for ind = 1:length(freqsReal)
  Pvals{ind} = Pappr(1i*freqsReal(ind));
end

epsgain = [10,50];
% epsgain = 13;
[ContrSys,epsgain] = PassiveRC(freqsReal,Pvals,epsgain,Sys);


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



CLSys = ConstrCLSys(Sys,ContrSys);
PlotEigs(full(CLSys.Ae))

stabmarg = CLStabMargin(CLSys)

xe0 = [x0;zeros(size(ContrSys.G1,1),1)];

Tend = 16;
tgrid = linspace(0,Tend,300);



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
hold off
cla
plot(tgrid,CLsim.error,'Linewidth',2);
set(gca,'xgrid','on','ygrid','on','tickdir','out','box','off')
title('Tracking  error $y(t)-y_{ref}(t)$','Interpreter','latex','Fontsize',16)
%set(gcf,'color',1/255*[252 247 255])


figure(3)
colormap jet
%PlotHeat2DSurf(x0,spgrid,[-1.4,1.4])
PlotLHMNCSurf(CLsim.xesol,spgrid,tgrid,[-9 9])
%PlotLHMNCSurf(CLsim.xesol(205:end,:),spgrid,tgrid,[-9 9])
set(gca,'ztick',-8:4:8);
%colormap jet

% figure(4)
% tt = linspace(0,16,500)
% plot(tt,yref(tt),'Color',1.1*[0 0.447 0.741],'Linewidth',3);
% set(gca,'xgrid','on','ygrid','on','tickdir','out','box','off')

