function CLSim = FreqEstimatorRC(Sys,x0,yref,wdist,tgrid,ContrConstructor,freqs_init,FreqEstMethod,T_step,T_freqs,T_end,freq_tol)
%
% A robust control scheme where the internal model is initialized with
% 'freqs_init', and the rest of the frequencies are estimated from the
% regulation error during the interval [0,T_freqs] (or until convergence of
% the frequencies is reached).
%
% freqs_init = the initial list of frequencies (to be used for 0<t<T_freqs
% ContrConstructor = A matlab for constructing the desired controller type
% based on the information on the list of frequencies (and possibly Pvals?)
% T_freqs = last limit after which the frequencies are added to IM (warning
% if convergence up to 'freq_tol' hasn't been reached
% T_step = Time-interval for checking the the convergence condition 
% T_end = end of the simulation
% freq_tol = tolerance for determining convergence of frequencies


freqs_init = 0;
FreqEstMethod = 'CarAst';

ContrSys = ContrConstructor(freqs_init);

dimX = size(Sys.A,1);
dimZ = size(ContrSys.G1,1);

dimXe = dimX + dimZ;


% The controller initial state is set to zero by default.
xe0 = [x0;zeros(dimZ,1)];

if isequal(FreqEstMethod,'FFT')
elseif isequal(FreqEstMethod,'CarAst')
    
    % Parameters of the frequency estimator
    % Multiple frequency case, CarAst08, Proposition 2:
    % number of frequencies to indentify
    q0 = 3;
    
    gamma1 = .005;
    gamma2 = 10;
    
    % parameter 'k': [k_1,...,k_{2*q0-1}] - Coefficients of a monic Hurwitz polynomial
    % of order 2*q0 (i.e., x^(2q0)+k_1*x^(2q0-1)+...+k_{2q0-1}
    kparam = poly(-2*ones(1,2*q0-1)); % Polynomial with all roots at s=-2
    kparam = kparam(2:end);
   
    % Construction of the matrices of the frequency estimator
    
    
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

    % FE_init = zeros(2*q0-1+q0,1);
    % Initialise the frequency parameters so that the corresponding frequencies
    % are nonzero and distinct.
    % Initial nonzero frequencies = [1,2,...,q0]
    est_freqs = 1:q0;
    
    
    % Compute the corresponding initial state for the frequency estimator (such
    % that xi_1(0)=0)
    theta0 = poly(-est_freqs.^2);
    theta0 = theta0(2:end).';
    FE_init = [zeros(2*q0-1,1);gamma1*theta0];
    
    % The full simulated dynamics consist of the closed-loop state xe
    % (dimXe) and the state xi of the frequency estimator.
    
    % The time-derivative of the frequency estimator part:
    CarAst_ODEfun = @(yaux,xi) [FE_Mainmat*xi(1:(2*q0-1))+(1/gamma1)*FE_B*yaux;-gamma2*Rfun(xi(1:(2*q0-1))).'*(FE_C.'*FE_C)*(Rfun(xi(1:(2*q0-1)))*(gamma2*yaux*(Rfun(xi(1:(2*q0-1))).'*FE_C.')+xi((2*q0):end))+FE_K*yaux+gamma1*Gmultmat*xi(1:(2*q0-1)))-gamma2*yaux*Rfun(FE_Mainmat*xi(1:(2*q0-1))+(1/gamma1)*FE_B*yaux).'*FE_C.'];

    % The time-derivative of the full simulated dynamics. The signal 'yaux'
    % for frequency estimate is the regulation error e(t)=y(t)-yref(t).
    CL_ODEfun = @(t,xee) [CLSys.Ae*xee(1:dimXe)+CLSys.Be*[wdist(t);yref(t)]; ... 
        CarAst_ODEfun(CLSys.Ce*xee(1:dimXe)+CLSys.De*[wdist(t);yref(t)],xee((dimXe+1):end))];
    
    xee0 = [xe0;FE_init];
    
    
elseif isequal(FreqEstMethod,'VedVed')
else
    error('Unknown frequency Estimation method.')
end    


% Simulate on [0,T_freq]
% Determine if the frequencies converged (up to freqs_tol)
% Construct a new controller with the frequency estimates

ContrSys_new = ContrConstructor(sort(unique([freqs_init,freqs_estimated])));


% Simulate on [T_freq,T_end]

CLSys_new = ContrCLSys(Sys,ContrSys_new);
CLSim_new = SimCLSys(CLSys_new,xe0,yref,wdist,tgrid_new);


% Compile the full results of the simulation


% IMPROVEMENT: THIS COULD BE AN ITERATIVE PROCESS
% - The frequency estimation could be continued, and new frequencies could
% be possibly added to the IM periodically after suitable intervals
% - Conversely, the "necessity" of the current frequencies could be
% evaluated during simulation (what would be a good criterion? The size of
% the controller state variables corresponding to these frequencies?)
% Mainly this process would be necessary if the frequency estimator
% identifies frequencies close to existing ones (in which case, these might
% be more accurate than earlier estimates).
% - The user could set, e.g., a maximum size of the internal model, or the
% desired limit for the tracking error after which the tuning process is
% considered to be "complete". Also max frequency could be good, to rule
% out estimation of noise. This could be helped with some filtering...

end


% % The closed-loop system ODE function for the 'CarAst'-method. 
% % The state xee = [xe,xi], where xe is the
% % state of the closed-loop system and xi is the state of the frequency
% % estimator 
% function dxidt = CarAst_ODEfun(yaux,xi,q0,FE_Mainmat,gamma1,gamma2,FE_C,FE_B,FE_K,Rfun,Gmultmat)
% 
% 
% % COMMENT OUT IF FREQS UPDATED AT INTERVALS
% % % frequency parameters
% % 
% % FE_out = (1/gamma1)*(gamma2*yaux*Rfun(xi(1:(2*q0-1))).'*FE_C.'+xi((2*q0):end));
% % % freqs = sort(sqrt(-roots([1,FE_out.'])));
% % freqs = [sqrt(abs(FE_out(1)+sqrt(abs(FE_out(1).^2-4*FE_out(2))))/2);sqrt(abs(FE_out(1)-sqrt(abs(FE_out(1).^2-4*FE_out(2))))/2)];
% 
% % % DEBUG Constant frequencies
% % freqs = [1,6];
% 
% 
% dxidt = [FE_Mainmat*xi(1:(2*q0-1))+(1/gamma1)*FE_B*yaux;-gamma2*Rfun(xi(1:(2*q0-1))).'*(FE_C.'*FE_C)*(Rfun(xi(1:(2*q0-1)))*(gamma2*yaux*(Rfun(xi(1:(2*q0-1))).'*FE_C.')+xi((2*q0):end))+FE_K*yaux+gamma1*Gmultmat*xi(1:(2*q0-1)))-gamma2*yaux*Rfun(FE_Mainmat*xi(1:(2*q0-1))+(1/gamma1)*FE_B*yaux).'*FE_C.'];
% 



end
