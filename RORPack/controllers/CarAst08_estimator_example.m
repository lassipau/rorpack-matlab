% Afshar-Paunonen 2020
% A 1D heat equation - output regulation with a time-dependent
% internal model


clear variables

global dimX dimY dimU dimUd dimZ0



% reference and disturbance signals
yref = @(t) 0.2*sin(0.5*t+0.5) + 0.4*sin(6*t+0.5);
wdist = @(t) [cos(1.5*t+0.5);sin(0.5*t+0.2);cos(6*t-0.4)];
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
Atilde= A + 0.5*eye(N);


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
% Multiple frequency case, CarAst08, Proposition 2:
% number of frequencies to indentify 
q0 = 3;

gamma1 = .005;
gamma2 = 10;

% parameter 'k': [k_1,...,k_{2*q0-1}] - Coefficients of a monic Hurwitz polynomial
% of order 2*q0 (i.e., x^(2q0)+k_1*x^(2q0-1)+...+k_{2q0-1}
kparam = poly(-2*ones(1,2*q0-1)); % Polynomial with all roots at s=-2
kparam = kparam(2:end);


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

figure(6)
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




%% The frequency estimator 

FE_Mainmat = flip(flip(compan([1,kparam]),1),2);
FE_C = [1,zeros(1,2*q0-2)];
FE_B = [zeros(2*q0-2,1);1];
FE_K = kparam.';

% Implement G(xi) as a multiplication Gmultmat*xi
Gmultmat1 = -kparam(:)*kparam(end:-1:1);
Gmultmat2 = zeros(2*q0-1);
for ind = 1:size(Gmultmat2,1)
    Gmultmat2(ind,:) = [zeros(1,ind),kparam(end:-1:(ind+1))];
end

Gmultmat = Gmultmat1+Gmultmat2;


Rmultmat = zeros(2*q0-1,q0,2*q0-1);
for indj = 1:size(Rmultmat,2)
    for indi = 1:(2*indj-1)
        % CASE i<2*j
        Rmultmat(indi,indj,:) = -[zeros(1,2*(q0-indj)),kparam((indi-1):-1:1),1,zeros(1,2*indj-1-indi)];
    end
    for indi = (2*indj):size(Rmultmat,1)
        % CASE i>=2*j
        Rmultmat(indi,indj,:) = [zeros(1,indi-2*indj),kparam((2*q0-1):-1:indi),zeros(1,2*indj-1)];
    end
    
end

% The function R(xi) from CarAst08
Rfun = @(xi) squeeze(pagemtimes(permute(Rmultmat,[1,3,2]),xi));


%%

Tmax = 20;
timestep = 1;
Nsteps = floor(Tmax/timestep);


xe0 = [x0;zeros(dimZ0+dimX,1)];

% FE_init = zeros(2*q0-1+q0,1);
% Initialise the frequency parameters so that the corresponding frequencies
% are nonzero and distinct.
% Initial nonzero frequencies = [1,2,...,q0]
est_freqs = 1:q0;
% est_freqs = [1,3,6];
% Compute the corresponding initial state for the frequency estimator (such
% that xi_1(0)=0)
theta0 = poly(-est_freqs.^2);
theta0 = theta0(2:end).';
FE_init = [zeros(2*q0-1,1);gamma1*theta0];
% FE_init = zeros(3*q0-1,1);

% Parameters for stabilization of the non-autonomous closed-loop system
Ae0 = @(freqs) [G1(freqs),G2(length(freqs))*C;zeros(dimX,dimZ0),A];
Be0 = @(freqs) [G2(length(freqs))*D;B];

% Parameters for the closed-loop stabilization
beta0 = .2;
Q = 1*eye(dimZ0+dimX);
R = eye(dimU);
% The minimal required stability margin for not updating the gain K
CL_stabmarg_par = 0.3;

% Tolerance parameter for overlapping frequencies and closeness to
% transmission zeros of (A,B,C,D).
eps_omega = 0.3;

%%
% Initial state of the combined closed-loop system and frequency estimator
xee0 = [xe0;FE_init];

ODE_opts = odeset('AbsTol',1e-12,'Reltol',1e-9);

Nt = 301; % Points on a single timestep
tt_step = linspace(0,timestep,Nt);

tt_full = zeros(1,Nsteps*(Nt-1)+1);
xe_full = zeros(dimX+dimZ0+dimX,Nsteps*(Nt-1)+1);
xi_full = zeros(3*q0-1,Nsteps*(Nt-1)+1);
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
            K1 = -lqr(G1(IM_freqs)+beta0*eye(dimZ0),B1j,10*eye(dimZ0),0.01*eye(dimU),zeros(dimZ0,dimU));

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

    sol = ode15s(@(t,xee) CL_odefun(t,xee,Ae,Be,Ceaux,Deaux,IM_freqs,K,dimX,dimZ0,wdist,yref,q0,FE_Mainmat,gamma1,gamma2,FE_C,FE_B,FE_K,Rfun,Gmultmat),[T,T+timestep],xee0,ODE_opts);
    
    xee_step = deval(sol,T+tt_step);
    % Range of indices (REMOVE REPEATING INDICES!!)
    if ind == 1
        indran = 1:Nt;
        tt_full(:,indran) = T+tt_step;
        xe_full(:,indran) = xee_step(1:(dimX+dimZ0+dimX),:);
        xi_full(:,indran) = xee_step((dimX+dimZ0+dimX+1):end,:);

        yy_full(:,indran) = Ce(K)*xee_step(1:(dimX+dimZ0+dimX),:);
        yaux_full(1,indran) = Ceaux*xee_step(1:(dimX+dimZ0+dimX),:)+Deaux*[wdist(tt_step);yref(tt_step)];
    else
        indran = (ind-1)*(Nt-1)+(2:Nt);
        tt_full(:,indran) = T+tt_step(2:end);
        xe_full(:,indran) = xee_step(1:(dimX+dimZ0+dimX),2:end);
        xi_full(:,indran) = xee_step((dimX+dimZ0+dimX+1):end,2:end);

        yy_full(:,indran) = Ce(K)*xee_step(1:(dimX+dimZ0+dimX),2:end);
        yaux_full(1,indran) = Ceaux*xee_step(1:(dimX+dimZ0+dimX),2:end)+Deaux*[wdist(T+tt_step(2:end));yref(T+tt_step(2:end))];

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

    
    % Compute the estimated frequencies at the time instances
    for ind2 = indran
        FE_out = (1/gamma1)*(gamma2*yaux_full(1,ind2)*Rfun(xi_full(1:(2*q0-1),ind2)).'*FE_C.'+xi_full((2*q0):end,ind2));
        % freqs_full(:,ind2) = sqrt(abs(FE_out(1)+[1;-1]*sqrt(abs(FE_out(1).^2-4*FE_out(2))))/2);
        freqs_full(:,ind2) = sort(sqrt(-roots([1,FE_out.'])));
    end
    
    
    xee0 = xee_step(:,end);
    
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
figure(3)
plot(tt_full,yaux_full,'Linewidth',2);
title('The residual output $y_{aux}(t)$.','Interpreter','Latex')
figure(4)
plot(tt_full,sqrt(sum(abs(yy_full-yref(tt_full)).^2,1)),'Linewidth',2);
% axis([tt_full(1),tt_full(end),0,0.2])
title('The norm of the regulation error $\|e(t)\|$.','Interpreter','Latex')

figure(5)
colormap jet
surf_indskip = 6;
Plot1DHeatSurf(xe_full(1:N,1:surf_indskip:end),spgrid,tt_full(1:surf_indskip:end),BCtype)

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
function dxeedt = CL_odefun(t,xee,Ae,Be,Ceaux,Deaux,IM_freqs,K,dimX,dimZ0,wdist,yref,q0,FE_Mainmat,gamma1,gamma2,FE_C,FE_B,FE_K,Rfun,Gmultmat)

xe = xee(1:(dimX+dimZ0+dimX));
xi = xee((dimX+dimZ0+dimX+1):end);

yaux = Ceaux*xe+Deaux*[wdist(t);yref(t)];

% COMMENT OUT IF FREQS UPDATED AT INTERVALS
% % frequency parameters
% 
% FE_out = (1/gamma1)*(gamma2*yaux*Rfun(xi(1:(2*q0-1))).'*FE_C.'+xi((2*q0):end));
% % freqs = sort(sqrt(-roots([1,FE_out.'])));
% freqs = [sqrt(abs(FE_out(1)+sqrt(abs(FE_out(1).^2-4*FE_out(2))))/2);sqrt(abs(FE_out(1)-sqrt(abs(FE_out(1).^2-4*FE_out(2))))/2)];

% % DEBUG Constant frequencies
% freqs = [1,6];

dxedt = Ae(IM_freqs,K)*xe + Be(IM_freqs)*[wdist(t);yref(t)];

dxidt = [FE_Mainmat*xi(1:(2*q0-1))+(1/gamma1)*FE_B*yaux;-gamma2*Rfun(xi(1:(2*q0-1))).'*(FE_C.'*FE_C)*(Rfun(xi(1:(2*q0-1)))*(gamma2*yaux*(Rfun(xi(1:(2*q0-1))).'*FE_C.')+xi((2*q0):end))+FE_K*yaux+gamma1*Gmultmat*xi(1:(2*q0-1)))-gamma2*yaux*Rfun(FE_Mainmat*xi(1:(2*q0-1))+(1/gamma1)*FE_B*yaux).'*FE_C.'];


dxeedt = [dxedt;dxidt];

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
Bdfun = @(xi) sin(3.5*xi);
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