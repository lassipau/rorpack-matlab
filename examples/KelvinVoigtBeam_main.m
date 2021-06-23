% Solve the 1D beam equation with Kelvin-Voigt damping on [-1,1] with clamped BCs
% using the Galerkin method with the Chebyshev function basis introduced in
% [Shen 1995]
% Copyright (C) 2020 by Lassi Paunonen (lassi.paunonen@tuni.fi)
% Licensed under GNU GPLv3 (see LICENSE.txt).


% Basis functions are \phi_k = T_k-2*(k+2)/(k+3)*T_{k+2}+(k+1)/(k+3)*T_{k+4} for k=0..N-4, where T_k are the
% Chebyshev polynomials of the first kind

% The inner products are defined as 
% <f,g>_w = \int_{-1}^1 f(x)g(x)(1-x^2)^{-1}dx

% Compile the approximation as a system M*v''(t) + a*F*v(t)+b*F*v'(t)+c*M*v'(t)= + Bu(t)
% Single output, B=b
% The formulas for the elements M_{kj} = (<\phi_j,\phi_k>_w) 
% and A_{kj} = (<\phi_j'',(w\phi_k)''>)

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
bd1 = @(r) (r+1).^2.*(1-r).^2;
dimUd = 1; % number of disturbance signals

% The initial profile (v0) and velocity (v0dot)
% v0fun = @(r) (r+1)^2.*(1-r).^3;
% v0fun = @(r) (r+1)^2.*(1-r).^2;
v0 = @(r) zeros(size(r));
v0dot = @(r) zeros(size(r));



% Size of the higher order approximation - used for simulation, representing the original PDE
Nhi = 70; 
% Size of the lower order approximation - used for controller design
Nlow = 40;



% Construct the system (the higher-dimensional Galerkin approximation)
[x0,Sys,Q_coeff] = ConstrEBKVbeam(E,I,d_KV,d_v,b1,b2,xi1,xi2,bd1,v0,v0dot,Nhi);



% Definition of the reference signal y_{ref}(t) and the disturbance signal
% w_{dist}(t), and the frequencies in the signal
yref1_oneper = @(t) 0.3+0.4*((t<=1).*(t-1/2)+(t>1).*(3/2-t));
yref = @(t) [yref1_oneper(rem(t,2));zeros(size(t))];
wdist = @(t) sin(pi*t)+0.4*cos(3*pi*t);
freqs = [pi*(0:10)];
q=10;

% yref = @(t) [yref1_oneper(rem(t,2));0.1*sin(2*t)];
% wdist = @(t) sin(pi*t)+0.4*cos(3*pi*t);
% freqs = [pi*(0:10),2];

% tt = linspace(0,10);
% plot(tt,yref(tt));



% % Parameters of the controller
% Galerkin approximation for the controller design

[~,Sys_Nlow,~] = ConstrEBKVbeam(E,I,d_KV,d_v,b1,b2,xi1,xi2,bd1,v0,v0dot,Nlow);

% Store the Galerkin approximation in "SysApprox".
SysApprox.AN = Sys_Nlow.A;
SysApprox.BN = Sys_Nlow.B;
SysApprox.CN = Sys_Nlow.C;
SysApprox.D = Sys_Nlow.D;

% Target stability margin of the closed-loop system (for LQR/LQG)
% CLstabmarg = 1;
alpha1 = 2;
alpha2 = 0.8;

Q0 = eye(IMdim(freqs,size(SysApprox.CN,1))); % Size = dimension of the IM 
Q1 = eye(size(SysApprox.AN,1)); % Size = dim(V_N)
Q2 = eye(size(SysApprox.AN,1)); % Size = dim(V_N)
R1 = eye(size(SysApprox.CN,1)); % Size = dim(Y)
R2 = eye(size(SysApprox.BN,2)); % Size = dim(U)

% Order of the reduced order observer in the controller
ROMorder = 4;


%%

% Parameters of the simulation
tspan = [0,16]; % time-interval of the simulation
% tspan = [0,50]; % time-interval of the simulation FOR LOW-GAIN

% Initial deflection profile (v_0) and the initial velocity (\dot{v}_0) 
% Initial state of the controller is by default zero


% Parameters of the visualisation of the results
% temporal grid for the output plot
tt_output = linspace(tspan(1),tspan(2),601);

% spatial and temporal grids for the state plot
spgrid_state = linspace(-1,1,81);
tt_state = linspace(tspan(1),tspan(2),201);

% spatial and temporal grids for the animation
spgrid_anim = linspace(-1,1,130);
tt_anim = linspace(tspan(1),tspan(2),401);
anim_pause = 0.02;



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


%% Construct the controller using a Galerkin approximation with size 'Nlo'


% Construct the Internal Model Based Reduced Order Controller
ContrSys = ConstrContrObsBasedROM(freqs,SysApprox,alpha1,alpha2,R1,R2,Q0,Q1,Q2,ROMorder);



% % For comparison: Construct the Low-Gain Internal Model Based Controller
% % Compute the values P(iw_k) using the Galerkin approximation of the 
% % control system
% % Adjust the set of frequencies, use only k\pi for k=0,...,5
% freqs = [pi*(0:5)];
% q=5;
% 
% Pvals = cell(length(freqs),1);
% for ind = 1:length(freqs)
%   Pvals{ind} = Sys_Nlo.C*((1i*freqs(ind)*eye(size(Sys_Nlo.A))-Sys_Nlo.A)\Sys_Nlo.B);
% end
% epsgain = 0.076;
% % epsgain = [0.01,0.2]; % The algorithm can optimize CL stability margin
% [ContrSys,epsgain] = ConstrContrLGReal(freqs,Pvals,epsgain,Sys_Nlo);


%% Closed-loop simulation


CLSys = ConstrCLSys(Sys,ContrSys);

stabmarg = CLStabMargin(CLSys)

figure(1)
PlotEigs(CLSys.Ae,[-20 .3 -150 150])

%%
xe0 = [x0;zeros(size(ContrSys.G1,1),1)];

Tend = 16;
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
plotControl(tgrid,CLsim,ContrSys,size(Sys.A,1),PrintFigureTitles)



%% Plot the state of the controlled beam equation
% Compute the values of the solution on the spatial grid tt_state
% by converting the coefficients (alpha_k) (given in the solution 'alphas') to
% the Chebyshev coeffients, and using them to define a Chebfun object at
% each timestep

xxvals_state = zeros(length(spgrid_state),length(tt_state));

xxe_state = deval(CLsim.solstruct,tt_state);
Cheb_coeffs = Q_coeff*xxe_state(1:(Nhi-1),:); % Position
% Cheb_coeffs = Q_coeff*xxe(Nhi:(2*Nhi-2),:); % Velocity
for ind = 1:length(tt_state)
  cfun = chebfun(Cheb_coeffs(:,ind),'coeffs');
  xxvals_state(:,ind) = cfun(spgrid_state);
end

figure(4)
surf(tt_state,spgrid_state,xxvals_state)
set(gca,'ydir','reverse')
ylabel('$\xi$','Interpreter','Latex','fontsize',18)
xlabel('$t$','Interpreter','Latex','fontsize',18)
if PrintFigureTitles, title('The deflection of the controlled beam','Interpreter','Latex','fontsize',16), end


%% Animate the solution of the controlled beam equation
xxvals_anim = zeros(length(spgrid_anim),length(tt_anim));

xxe_anim = deval(CLsim.solstruct,tt_anim);
Cheb_coeffs = Q_coeff*xxe_anim(1:(Nhi-1),:); % Position
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