% Afshar-Paunonen 2020
% A 1D heat equation - output regulation with a time-dependent
% internal model
% Frequency estimation using the finite time adaptive estimator in
% VedVed20arxiv


clear variables

global dimX dimY dimU dimUd dimZ0


true_freqs = sort([0.5;1.5;6]);


% reference and disturbance signals
yref = @(t) 0.2*sin(0.5*t+0.5) + 0.4*sin(6*t+0.5);
wdist = @(t) [cos(1.5*t+0.5);sin(0.5*t+0.2);cos(1.5*t-0.4)];
% yaux = @(t) sin(2*t) + sin(5*t);
% yaux = @(t) sin(2*t); %+ 0.4*exp(-0.4*t);
% yaux = @(t) sin(2*t) + sin(5*t) + sin(7*t);

% Initial temperature profile
x0fun = @(xi) zeros(size(xi));


% Define the model
N = 100; % Dimension of the FD approximation
[x0,Sys,spgrid,BCtype] = Constr1DHeatModel(1,x0fun,N);

A = Sys.A;
B = Sys.B;
Bd = Sys.Bd;
C = Sys.C;
D = Sys.D;

% Unmodelled reaction term as a perturbation
Atilde= A + 0.001*eye(N);


% % The system
% A = [-1,1;0,2];
% B = [1,1;1,2];
% Bd = [0.5;0];
% C = [2,1;1,2];
% D = [0.05,0;0,0.1];

% max(eig(A-100*(C'*C)))
 
% Check that the system does not have transmission zeros on iR
TransZeros = tzero(ss(A,B,C,D));
TZ_ind = find(abs(real(TransZeros))<1e-6);
TZ_imag = abs(imag(TransZeros(TZ_ind)));

% Test purposes, consider two "transmission zero" points
% TZ_imag = [2.5,4].';



%%

% The parameters of the frequency estimator
% Multiple frequency case, VedVed20arxiv:
% In this implementation the frequency estimator is not run "online", but
% instead it is run one interval at a time at the update times. The online
% implementation would be a bit more challenging since the estimator uses
% delayed values of the y_aux(t) signal.

% number of frequencies to indentify 
q0 = 3;

eps_par = 100;
gamma_par = 10*ones(q0,1);
h = 0.1;
d = 0.2;

delay_values = h*(0:(2*q0)).'*ones(1,q0)+ones(2*q0+1,1)*d*(1:q0);


% Coefficients C_n^i = nchoosek(n,n-i)
Cni = @(n_val,i_val) nchoosek(n_val,n_val-i_val);

% Matrix 'B' in VedVed20arxiv 
FE_B = zeros(q0,2*q0-1);
for ind = 1:q0
    Cvec_tmp = zeros(2,q0-ind+1);
    for ind2 = 1:(q0-ind+1)
        Cvec_tmp(:,ind2)=[Cni(q0-ind,q0-ind-(ind2-1));0];
    end
    Cvec_flat = Cvec_tmp(:).';
    Cvec_flat = Cvec_flat(1:(end-1));
    FE_B(ind,:) = 2^ind*[zeros(1,ind-1),Cvec_flat,zeros(1,ind-1)];
end

% Construct a matrix FE_Phi_mat for computing
% Phi_f(t)=(FE_Phi_mat*yaux_vals(t)).'
% Defined as an extended FE_B, to be applied directly to yaux_vals(t) 
% (leaves out first and last row of yaux_vals(t))
FE_Phi_mat = [zeros(q0,1),FE_B,zeros(q0,1)];

% Construct a matrix FE_Psi_mat for computing
% Psi_f(t)=FE_Psi_mat*yaux_vals(t).
% Defined as [Cni(q0,q0),0,Cni(q0,q0-1),0,...,0,Cni(q0,0)]
FE_Psi_mat = zeros(2,q0+1);
for ind = 1:(q0+1)
    FE_Psi_mat(1,ind) = Cni(q0,q0-ind+1);
end
FE_Psi_mat = FE_Psi_mat(1:(end-1));


% Stabilisation of the observer

dimX = size(A,2);
dimY = size(C,1);
dimU = size(B,2);
dimUd = size(Bd,2);

% The number of estimated frequencies determines the size of the internal
% model
dimZ0 = (2*q0+1)*dimY;

% Choose L and K21 so that A+L*C and A+B*K21 are exponentially stable
% L = -place(A',C',[-1,-2])';
L = -4*ones(N,1);
K21 = 2*ones(1,N)/(N-1);

% figure(6)
% PlotEigs(full(A+L*C),[-20 1 -.3 .3])
% PlotEigs(full(A+B*K21),[-20 1 -.3 .3])



%%





% Internal model (TO IMPROVE!)
% G1 = @(freqs) blkdiag(0,freqs(1)*[0,1;-1,0],freqs(2)*[0,1;-1,0],freqs(3)*[0,1;-1,0]);
G1 = @(freqs) blkdiag(zeros(dimY),kron(diag(freqs),[zeros(dimY),eye(dimY);-eye(dimY),zeros(dimY)]));
% G2 = [1;1;0;1;0;1;0];
G2 = @(Nfreqs) [eye(dimY);kron(ones(Nfreqs,1),[eye(dimY);zeros(dimY)])];


G1full = @(freqs,K) [G1(freqs),zeros(dimZ0,dimX);zeros(dimX,dimZ0),A+L*C] + [zeros(dimZ0,dimZ0+dimX);(B+L*D)*K];
G2full = @(Nfreqs) [G2(Nfreqs);-L];


Ae = @(freqs,K) [Atilde,B*K;G2full(length(freqs))*C,G1full(freqs,K)+G2full(length(freqs))*D*K];
Be = @(freqs) [Bd,zeros(dimX,dimY);zeros(dimZ0+dimX,dimUd),-G2full(length(freqs))];

Ce = @(K) [C,D*K];
De = -eye(dimY);

% y_{aux}(t) as an output of the closed-loop system:
% y_{aux}(t) = [C,0,-C]*xe+[0,I]*[wdist(t),yref(t)]
% y_{aux}(t) is always scalar-valued, for multiple outputs, define y_aux(t)
% as a random linear combination of the original components.
if dimY > 1
    rgen = rng(1,'twister');
    % rand('seed',0);
    yaux_randvec = rand([1,dimY]);
else
    yaux_randvec=1;
end

Ceaux = yaux_randvec*[C,zeros(dimY,dimZ0),-C];
Deaux = yaux_randvec*[zeros(dimY,dimUd),-eye(dimY)];



%%

Tmax =20;
timestep = 1;
Nsteps = floor(Tmax/timestep);


xe0 = [x0;zeros(dimZ0+dimX,1)];



% Initialise the frequency estimator so that the corresponding frequencies
% are nonzero and distinct.
% Initial nonzero frequencies = [1,2,...,q0]
est_freqs = 1:q0;

ci_init = cos(h*est_freqs);
syms x 
init_pol = x-ci_init(1);
for ind_init = 2:q0
    init_pol = init_pol*(x-ci_init(ind_init));
end
FE_init = double(coeffs(init_pol,'All'));
FE_init = -FE_init(2:end).';



% Parameters for stabilization of the non-autonomous closed-loop system
Ae0 = @(freqs) [G1(freqs),G2(length(freqs))*C;zeros(dimX,dimZ0),A];
Be0 = @(freqs) [G2(length(freqs))*D;B];

% Parameters for the closed-loop stabilization
beta0 = .2;
Q1 = 10*eye(dimZ0);
R1 = 0.01*eye(dimU);

% The minimal required stability margin for not updating the gain K
CL_stabmarg_par = 0.3;

% Tolerance parameter for overlapping frequencies and closeness to
% transmission zeros of (A,B,C,D).
eps_omega = 0.3;

%%

ODE_opts = odeset('AbsTol',1e-12,'Reltol',1e-9);

Nt = 101; % Points on a single timestep
tt_step = linspace(0,timestep,Nt);

tt_full = zeros(1,Nsteps*(Nt-1)+1);
xe_full = zeros(dimX+dimZ0+dimX,Nsteps*(Nt-1)+1);
freqs_full = zeros(q0,Nsteps*(Nt-1)+1);
yy_full = zeros(dimY,Nsteps*(Nt-1)+1);

% y_{aux}(t) is always scalar-valued
yaux_full = zeros(1,Nsteps*(Nt-1)+1);

% % DEBUG Constant frequencies
% freqs = [1,6];

tic
for ind = 1:Nsteps
    
    
    % Current time
    T = (ind-1)*timestep;

    % Test whether or not the frequencies in the internal model should be
    % updated:
    
    % Variable 'true' if there are no overlapping frequencies (within the
    % tolerance 'eps_omega')
    freq_cond_overlap = (length(real(est_freqs))==length(uniquetol(real(est_freqs),eps_omega,'DataScale',1)));
    
    % Variable 'true' if none of the estimated frequencies are close to the
    % imaginary transmission zeros of (A,B,C,D) (within the tolerance 
    % 'eps_omega')
    freq_cond_TZ = ~any(ismembertol([real(est_freqs(:));-real(est_freqs(:))],TZ_imag,eps_omega,'DataScale',1));

    % Variable 'true' if all frequencies are real and positive
    freq_cond_real = ~any(imag(est_freqs)) && all(est_freqs>0);
    
    % If all the conditions are satisfied, update the frequencies in the
    % internal model. 
    if freq_cond_overlap && freq_cond_TZ && freq_cond_real
        IM_freqs = est_freqs;
    
        % Compute a new stabilizing gain if the old K does not stabilize the
        % pair (A_{e0},B_{e0}) with the updated frequencies (with a stability
        % margin of at least "CL_stabmarg_par"). Otherwise do not update K.
        if ind==1 || max(real(eig(full(Ae0(IM_freqs)+Be0(IM_freqs)*K))))>-CL_stabmarg_par
            
            % Compute the stabilizing gain K for the current timestep
            % Pi = are(Ae0(IM_freqs)+beta0*eye(dimZ0+dimX),Be0(IM_freqs)*(R\Be0(IM_freqs)'),Q);
            % K = -R\(Be0(IM_freqs)'*Pi);
            
            Hj = sylvester(G1(IM_freqs),-(A+B*K21),G2(length(IM_freqs))*(C+D*K21));
            B1j = Hj*B+G2(length(IM_freqs))*D;

            % Stabilization of the internal model, choose K1 so that G1+B1*K1 is
            % exponentially stable
            K1 = -lqr(G1(IM_freqs)+beta0*eye(dimZ0),B1j,Q1,R1,zeros(dimZ0,dimU));

            K = [K1,K21+K1*Hj];

        elseif ind~=1
            fprintf(['Parameter K not updated at time t = ' num2str(T) '\nClosed-loop stability margin = ' num2str(-max(real(eig(full(Ae0(IM_freqs)+Be0(IM_freqs)*K))))) '\n'])
            
        end
    elseif ind==1
        error('The initial frequency estimates are invalid');
    else
        fprintf(['Frequencies not updated at time t = ' num2str(T) '\nViolated conditions:\n'])
        if ~freq_cond_overlap, fprintf('Overlapping frequencies\n'),end
        if ~freq_cond_TZ, fprintf('Frequencies close to transmission zeros\n'),end
        if ~freq_cond_real, fprintf('Negative or complex frequency estimates\n'),end
        fprintf('\n')
    end

    sol = ode15s(@(t,xee) CL_odefun(t,xee,Ae,Be,Ceaux,Deaux,IM_freqs,K,dimX,dimZ0,wdist,yref),[T,T+timestep],xe0,ODE_opts);
    
    xe_step = deval(sol,T+tt_step);
    % Range of indices (REMOVE REPEATING INDICES!!)
    if ind == 1
        indran = 1:Nt;
        tt_full(:,indran) = T+tt_step;
        xe_full(:,indran) = xe_step(1:(dimX+dimZ0+dimX),:);

        yy_full(:,indran) = Ce(K)*xe_step(1:(dimX+dimZ0+dimX),:);
        yaux_full(1,indran) = Ceaux*xe_step(1:(dimX+dimZ0+dimX),:)+Deaux*[wdist(tt_step);yref(tt_step)];
    else
        indran = (ind-1)*(Nt-1)+(2:Nt);
        tt_full(:,indran) = T+tt_step(2:end);
        xe_full(:,indran) = xe_step(1:(dimX+dimZ0+dimX),2:end);

        yy_full(:,indran) = Ce(K)*xe_step(1:(dimX+dimZ0+dimX),2:end);
        yaux_full(1,indran) = Ceaux*xe_step(1:(dimX+dimZ0+dimX),2:end)+Deaux*[wdist(T+tt_step(2:end));yref(T+tt_step(2:end))];

    end
 
    
%     % Frequency estimation with FFT:
%     % First: Based on full history of y_{aux}(t)
% 
%     tau = tt_full(2)-tt_full(1);
%     N_freq = floor(indran(end)/2)*2; % make sure that the index is even
%     yaux_freq = yaux_full(:,1:N_freq);
%     yaux_fft = fft(yaux_freq.',N_freq).';
%     yaux_power = sum(abs(yaux_fft).^2,1)/N_freq;
%     freq_axis = 2*pi/(tau*N_freq)*(0:(N_freq/2-1));
% %     plot(freq_axis,yaux_power(1:(N_freq/2)))
%      xlim([0,12])
% %  
% %     ylim([0,15])
% 
% 
%     [freq_strengths,all_freq_inds]=findpeaks(yaux_power(1:(N_freq/2)),'SortStr','descend');
%     all_est_freqs = freq_axis(all_freq_inds);
% 
%     fft_freqs = sort(all_est_freqs(1:min(q0,length(all_est_freqs))));

    % VedVed20arxiv estimator
    % The delayed values of y_aux(t) are obtained from yaux_full with
    % linear interpolation. The values for t<0 are defined to be zero.
    
    yaux_vals = @(t) interp1([-q0*(2*h+d),-.001,tt_full(1:indran(end))],[0,0,yaux_full(1:indran(end))],t-delay_values);
    
    Phi_f = @(t) (FE_Phi_mat*yaux_vals(t)).';
    Psi_f = @(t) (FE_Psi_mat*yaux_vals(t)).';
    
    Psi = @(t) adjoint(eps_par*Phi_f(t))*eps_par*Psi_f(t);
    
    Delta = @(t) det(eps_par*Phi_f(t));
    FE_ode = @(t,theta_hat) Delta(t)*diag(gamma_par)*(Psi(t)-Delta(t)*theta_hat);

    
    opts = odeset('AbsTol',1e-9,'Reltol',1e-6);
    opts = odeset();
    FE_sol = ode15s(FE_ode,[T,T+timestep],FE_init,opts);

    theta_vals = deval(FE_sol,T+tt_step);

    FE_init = FE_sol.y(:,end);
    
    FE_W_ints = zeros(1,length(tt_step));
    for ind_FE = 2:length(tt_step)
        FE_W_ints(ind_FE) = FE_W_ints(ind_FE-1)+(tt_step(2)-tt_step(1))*(Delta(T+tt_step(ind_FE))^2+Delta(T+tt_step(ind_FE-1))^2)/2;
    end
    FE_W_vals = exp(-FE_W_ints);
    FE_W_gamma_vals = FE_W_vals.^gamma_par;
    
    % The first column of theta_ft_vals contains a division by zero,
    % insert the theoretical value for this column directly.
    theta_ft_vals = [theta_vals(:,1),(theta_vals(:,2:end)-diag(theta_vals(:,1))*FE_W_gamma_vals(:,2:end))./(1-FE_W_gamma_vals(:,2:end))];
    % UNCOMMENT TO USE \omega^{grad}
    % theta_ft_vals=theta_vals;

    omega_vals = zeros(q0,length(tt_step));
    for ind_FE = 1:length(tt_step)
        omega_vals(:,ind_FE) = 1/h*sort(acos(roots([1,-theta_vals(:,ind_FE).'])));
    end
    if ind == 1
        freqs_full(:,indran) = omega_vals;
    else
        freqs_full(:,indran) = omega_vals(:,2:end);
    end

    
    
    xe0 = xe_step(:,end);
    
    est_freqs = freqs_full(:,indran(end));


    
end
toc
%%

figure(1)
hold off
cla
hold on
plot(tt_full,yref(tt_full),'color',0.5*[1,1,1],'Linewidth',2)
plot(tt_full,yy_full,'Linewidth',2)
title('Regulated output $y(t)$.','Interpreter','Latex')
figure(2)
plot(tt_full,freqs_full,'Linewidth',2);
title('Estimated frequencies.','Interpreter','Latex')
% plot(tt_full,sqrt(sum(abs(freqs_full-true_freqs*ones(1,length(tt_full))).^2,1)),'Linewidth',2);
% title('Frequency estimation error (Euclidean norm).','Interpreter','Latex')
% plot(tt_full,max(abs(freqs_full-true_freqs*ones(1,length(tt_full))),[],1),'Linewidth',2);
% title('Frequency estimation error (Euclidean norm).','Interpreter','Latex')
figure(3)
plot(tt_full,yaux_full,'Linewidth',2);
title('The residual output $y_{aux}(t)$.','Interpreter','Latex')
figure(4)
plot(tt_full,sqrt(sum(abs(yy_full-yref(tt_full)).^2,1)),'Linewidth',2);
% axis([tt_full(1),tt_full(end),0,0.2])
title('The norm of the regulation error $\|e(t)\|$.','Interpreter','Latex')

figure(5)
% colormap jet
surf_t_indskip = 6;
surf_xi_indskip = 5;
Plot1DHeatSurf(xe_full(1:surf_xi_indskip:N,1:surf_t_indskip:end),spgrid(1:surf_xi_indskip:N),tt_full(1:surf_t_indskip:end),BCtype)


figure(6)
subplot(2,1,1)
plot(tt_full,freqs_full,'Linewidth',2);
title('Estimated frequencies $\hat{\omega}_k(t)$.','Interpreter','Latex')
% xlabel('t','Interpreter','Latex')
subplot(2,1,2)
plot(tt_full,sqrt(sum(abs(yy_full-yref(tt_full)).^2,1)),'Linewidth',2);
% axis([tt_full(1),tt_full(end),0,0.2])
title('The norm of the regulation error $\|e(t)\|$.','Interpreter','Latex')
% xlabel('t','Interpreter','Latex')


%%


% The closed-loop system ODE function, state xee = [xe,xi], where xe is the
% state of the closed-loop system and xi is the state of the frequency
% estimator
function dxedt = CL_odefun(t,xee,Ae,Be,Ceaux,Deaux,IM_freqs,K,dimX,dimZ0,wdist,yref)

xe = xee(1:(dimX+dimZ0+dimX));

% yaux = Ceaux*xe+Deaux*[wdist(t);yref(t)];


dxedt = Ae(IM_freqs,K)*xe + Be(IM_freqs)*[wdist(t);yref(t)];

end


% Construct the parameters of the heat equation model
function [x0,Sys,spgrid,BCtype] = Constr1DHeatModel(cval,x0fun,N)

spgrid = linspace(0,1,N);
h = 1/(N-1);

ee = ones(N,1);

A = cval*1/h^2*full(spdiags([ee -2*ee ee],-1:1,N,N));
A(1,2) = cval*2/h^2;
A(N,N-1) = cval*2/h^2;


% Neumann boundary input at x=0, sign is based on the _outwards_ normal derivative
B = [-2/h;zeros(N-1,1)]; 
Bd1 = [-2/h;zeros(N-1,1)]; % First Neumann boundary disturbance at x=0


Bd2 = [zeros(N-1,1);2/h]; % Second Neumann boundary disturbance at x=1
C = [zeros(1,N-1), 1]; % Measured temperature at x=1

% Case 1 has Neumann boundary conditions at both x=0 and x=1
BCtype = 'NN';

% Third distributed disturbance with (an unknown) profile
Bdfun = @(xi) cos(3*xi);
Bd3 = Bdfun(spgrid).';



x0 = x0fun(spgrid).';


Sys.A = A;
Sys.B = B;
Sys.Bd = [Bd1,Bd2,Bd3];
Sys.C = C;
Sys.D = 0;
Sys.Dd = 0;
end

function Plot1DHeatSurf(state,spgrid,tgrid,BCtype,zlims)
% function Plot1DHeatSurf(state,spgrid,tgrid,zlims,BCtype)
%
% Plot the solution (temperature) of the controlled 1D heat equation
% state = state of the heat system at times 'tgrid'. One column for each
% time instance
% spgrid = spatial grid
% tgrid = grid for time
% BCtype = 'NN' for Neumann-Neumann, 'ND' for Neumann-Dirichlet etc
% zlims = limits for the z-axis (optional)

if isequal(BCtype(1),'D')
  state = [zeros(1,length(tgrid));state];
elseif ~isequal(BCtype(1),'N')
  error('Unknown boundary condition type')
end

if isequal(BCtype(2),'D')
  state = [state;zeros(1,length(tgrid))];
elseif ~isequal(BCtype(2),'N')
  error('Unknown boundary condition type')
end

if max(max(abs(imag(state)))) > 1e-8
  warning('Solution may contain imaginary parts that are ignored in the plot.')
end

state = real(state);

if nargin <= 5
  zlims(1) = min(min(min(state)));
  zlims(2) = max(max(max(state)));
end
axlims = [tgrid(1) tgrid(end) spgrid(1) spgrid(end) zlims];
    
surf(tgrid,spgrid,state);
axis(axlims)
caxis(zlims)
set(gca, 'YDir', 'reverse')

xlabel('$t$','Interpreter','latex','Fontsize',20)
ylabel('$\xi$','Interpreter','latex','Fontsize',20)

end

function PlotEigs(A,axlim)
% function PlotEigs(A,axlims)
%
% Plots the eigenvalues of A
% If 'axlim' is not given, limits determined from the spectrum. If some
% components of 'axlim' are set as NaN, then also those limits will also be 
% determined from the spectrum.

Aspec = eig(full(A));

if nargin == 1
  relim = [min(real(Aspec)) max(real(Aspec))];
  imlim = [min(imag(Aspec)) max(imag(Aspec))];
  axlim = [relim imlim];
else
  if isnan(axlim(1))
    axlim(1) = min(real(Aspec));
  end
  if isnan(axlim(2))
    axlim(2) = max(real(Aspec));
  end
  if isnan(axlim(3))
    axlim(3) = min(imag(Aspec));
  end
  if isnan(axlim(4))
    axlim(4) = max(imag(Aspec));
  end
end
  

hold off
cla
hold on
plot(real(Aspec),imag(Aspec),'r.','Markersize',15)
% set the limits of the plot
axis(axlim)

% plot the axes
plot(axlim(1:2),[0 0],'k',[0 0],axlim(3:4),'k','Linewidth',1)


maxreal = num2str(max(real(Aspec)));
maxrealLF = num2str(max(real(Aspec(abs(imag(Aspec))<100))));
title(['Largest real part = $' maxreal '$ (Low freq = $ ' maxrealLF '$)' ],'Interpreter','Latex')
end