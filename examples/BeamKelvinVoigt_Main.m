% Robust internal model based controller design for a 1D Euler-Bernoulli 
% beam model with Kelvin-Voigt damping on [-1,1]. The beam has clamped
% boundary conditions, 2 distributed inputs and 2 pointwise measurements.
% The controller design is completed using a "reduced order internal model
% based controller" from the reference Paunonen & Phan IEEE TAC 2020. In
% RORPack the controller design is implemented in the routine
% 'ObserverBasedROMRC'.
% The system is approximated using the Galerkin method based on a nonlocal 
% basis based on Chebyshev polynomials (introduced in [Shen 1995]).
% 
% The simulation example is part of the conference paper "A Reduced 
% Order Controller for Output Tracking of a Kelvin–Voigt Beam" by Paunonen
% and Phan at MTNS 2020 (available at https://arxiv.org/abs/2104.08025)



% Physical parameters of the system
E = 10;
% E = 2.1e11;
I = 1;
% I = 1.167e-10;
% Damping coefficients d_KV (Kelvin-Voigt damping) and d_v (viscous
% damping)
% d_KV = .05;
d_KV = 0.01;
% d_v = 0.5;
% d_v = 0.001;
d_v = 0.4;

% Input profile functions
b1 = @(xi) 1/3*(xi+1).^2.*(1-xi).^6;
b2 = @(xi) 1/3*(xi+1).^6.*(1-xi).^2;

% Locations of the pointwise observations y(t)=[v(\xi_1,t);v(\xi_2,t)]^T
xi1 = -0.6;
xi2 = 0.3;

% Disturbance input profile function
bd1 = @(xi) (xi+1).^2.*(1-xi).^2;
dimUd = 1; % number of disturbance signals


% Initial deflection profile (v_0) and the initial velocity (\dot{v}_0) 
% v0fun = @(r) (r+1)^2.*(1-r).^3;
% v0fun = @(r) (r+1)^2.*(1-r).^2;
v0 = @(r) zeros(size(r));
v0dot = @(r) zeros(size(r));



% Order of the numerical approximation of the beam model. (This is the
% numerical approximation of the plant, which is used for simulation. The
% controller design uses a separate lower order Galerkin approximation.)
N = 70; 

% Construct the approximation of the PDE model 
[x0,Sys,Q_coeff] = ConstrBeamKelvinVoigt(E,I,d_KV,d_v,b1,b2,xi1,xi2,bd1,v0,v0dot,N);



% Definition of the reference signal y_{ref}(t) and the disturbance signal
% w_{dist}(t), and the frequencies in the signal
yref1_oneper = @(t) 0.3+0.4*((t<=1).*(t-1/2)+(t>1).*(3/2-t));
yref = @(t) [yref1_oneper(rem(t,2));zeros(size(t))];
wdist = @(t) sin(pi*t)+0.4*cos(3*pi*t);
freqsReal = pi*(0:10);

% yref = @(t) [yref1_oneper(rem(t,2));0.1*sin(2*t)];
% wdist = @(t) sin(pi*t)+0.4*cos(3*pi*t);
% freqs = [pi*(0:10),2];

% tt = linspace(0,10);
% plot(tt,yref(tt));

% %% Plot the norm of the inverse \|P(is)^{-1}\| of the transfer function on iR 
% % (to study the locations of possible transmission zeros)
% 
% ss = linspace(0,20);
% Ptransfun = @(s) C*((s*eye(2*(N-1))-A)\B);
% Pinvnorms = zeros(size(ss));
% for ind = 1:length(ss)
%   tmpval=svd(Ptransfun(1i*ss(ind)));
% % tmpval=real(Ptransfun(1i*ss(ind)));
%   Pinvnorms(ind) = 1/tmpval(2);
% end
% plot(ss,Pinvnorms)

% Check the consistency of the system definition
Sys = SysConsistent(Sys,yref,wdist,freqsReal);


%% Construct the reduced order controller 

% Parameters of the controller
% 
% Construct the Galerkin approximation for the controller design.
% The order 'Nlow' of the Galerkin approximation used in the controller 
% design lower than the approximation of the plant in 'Sys'.
Nlow = 40;
[~,Sys_Nlow,~] = ConstrBeamKelvinVoigt(E,I,d_KV,d_v,b1,b2,xi1,xi2,bd1,v0,v0dot,Nlow);

% Store the Galerkin approximation in "SysApprox".
SysApprox.AN = Sys_Nlow.A;
SysApprox.BN = Sys_Nlow.B;
SysApprox.CN = Sys_Nlow.C;
SysApprox.D = Sys_Nlow.D;

% Target stability margin of the closed-loop system in the controller 
% design and other design parameters (for LQR/LQG design)
alpha1 = 2;
alpha2 = 0.8;

Q0 = eye(IMdim(freqsReal,size(SysApprox.CN,1))); % Size = dimension of the IM 
Q1 = eye(size(SysApprox.AN,1)); % Size = dim(V_N)
Q2 = eye(size(SysApprox.AN,1)); % Size = dim(V_N)
R1 = eye(size(SysApprox.CN,1)); % Size = dim(Y)
R2 = eye(size(SysApprox.BN,2)); % Size = dim(U)

% Order of the reduced order observer in the controller
ROMorder = 4;

% Construct the Internal Model Based Reduced Order Controller
ContrSys = ObserverBasedROMRC(freqsReal,SysApprox,alpha1,alpha2,R1,R2,Q0,Q1,Q2,ROMorder);


% % For comparison: Construct the Low-Gain Internal Model Based Controller
% % Compute the values P(iw_k) using the Galerkin approximation of the 
% % control system
% % Adjust the set of frequencies, use only k\pi for k=0,...,5
% freqs = [pi*(0:5)];
% q=5;
% 
% Compute the transfer function values in 'Pvals'. Since the controller
% feedthrough will be set to zero (Dc=0), the elements in Pvals are the
% values P(1i*w_k) of the transfer function of (A,B,C,D) where w_k are the
% frequencies in 'freqsReal'. Here we also use a lower order approximation
% 'Sys_Nlow' for the controller construction.
%
% Nlow = 40;
% [~,Sys_Nlow,~] = ConstrBeamKelvinVoigt(E,I,d_KV,d_v,b1,b2,xi1,xi2,bd1,v0,v0dot,Nlow);
%
% Pvals = cell(length(freqs),1);
% for ind = 1:length(freqs)
%   Pvals{ind} = Sys_Nlo.C*((1i*freqs(ind)*eye(size(Sys_Nlo.A))-Sys_Nlo.A)\Sys_Nlo.B);
% end
% epsgain = 0.076;
% % epsgain = [0.01,0.2]; % The algorithm can optimize CL stability margin
% [ContrSys,epsgain] = LowGainRC(freqs,Pvals,epsgain,Sys_Nlo);


%% Closed-loop construction and simulation

% Construct the closed-loop system
CLSys = ConstrCLSys(Sys,ContrSys);

% Print an approximate stability margin of the closed-loop system
stabmarg = CLStabMargin(CLSys)

% Define the initial state of the closed-loop system
% (the controller has zero initial state by default).
xe0 = [x0;zeros(size(ContrSys.G1,1),1)];

% Set the simulation length and define the plotting grid
Tend = 14;
% Tend = 50; % A longer time-interval for the LOW-GAIN CONTROLLER
tgrid = linspace(0,Tend,300);

% Simulate the closed-loop system
CLsim = SimCLSys(CLSys,xe0,yref,wdist,tgrid,[]);

%% Visualization

% Plot the eigenvalues of the (approximate) closed-loop system
figure(1)
PlotEigs(CLSys.Ae,[-20 .3 -150 150])

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
set(gcf,'color','white')

%% Plot the state of the controlled beam equation
% Compute the values of the solution on the spatial grid tt_state
% by converting the coefficients (alpha_k) (given in the solution 'alphas') to
% the Chebyshev coeffients, and using them to define a Chebfun object at
% each timestep

% spatial and temporal grids for the state plot
spgrid_state = linspace(-1,1,81);
tt_state = linspace(0,Tend,201);

xxvals_state = zeros(length(spgrid_state),length(tt_state));

xxe_state = deval(CLsim.solstruct,tt_state);
Cheb_coeffs = Q_coeff*xxe_state(1:(N-1),:); % Position
% Cheb_coeffs = Q_coeff*xxe(Nhi:(2*Nhi-2),:); % Velocity
for ind = 1:length(tt_state)
  cfun = chebfun(Cheb_coeffs(:,ind),'coeffs');
  xxvals_state(:,ind) = cfun(spgrid_state);
end


figure(4)
surf(tt_state,spgrid_state,xxvals_state,'Linewidth',0.2)
set(gca,'ydir','reverse')
ylabel('$\xi$','Interpreter','Latex','fontsize',18)
xlabel('$t$','Interpreter','Latex','fontsize',18)
if PrintFigureTitles, title('The deflection of the controlled beam','Interpreter','Latex','fontsize',16), end


%% Animate the solution of the controlled beam equation

% spatial and temporal grids for the animation
spgrid_anim = linspace(-1,1,130);
tt_anim = linspace(0,Tend,401);
anim_pause = 0.02;

xxvals_anim = zeros(length(spgrid_anim),length(tt_anim));

xxe_anim = deval(CLsim.solstruct,tt_anim);
Cheb_coeffs = Q_coeff*xxe_anim(1:(N-1),:); % Position
% Cheb_coeffs = Q_coeff*xxe_animation(Nhi:(2*Nhi-2),:); % Velocity
for ind = 1:length(tt_anim)
  cfun = chebfun(Cheb_coeffs(:,ind),'coeffs');
  xxvals_anim(:,ind) = cfun(spgrid_anim);
end


figure(5)
% No movie recording
[~,zlims] = Animate1DKVBeam(xxvals_anim,spgrid_anim,tt_anim,anim_pause,0);

% Movie recording
% [MovAnim,zlims] = Animate1D(xxvals_anim,spgrid_anim,tt_anim,anim_pause,1);
