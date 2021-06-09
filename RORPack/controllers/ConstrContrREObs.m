function [ContrSys,K_S] = ConstrContrREObs(freqs,Sys,K_S,L1,gammafun)
% ContrSys = ConstrContrObsBasedReal(freqs,Sys)
%
% Construct a robust controller based on solving the regulator equations for SISO systems
% The frequencies are assumed to be conjugate pairs, and the internal 
% model is in real form
% freqs = Frequencies to be included in the controller, only real nonnegative
% frequencies, if zero frequency is included, it's the first element in the
% vector. The control system is assumed to be real (i.e.,
% P(conj(s))=conj(P(s))), and P(iw_k) are invertible at the frequencies of
% the reference and disturbance signals.
% Pvals = [cell array] Values (or approximations of them) of the values of the transfer
% function of the system on the frequencies 'freqs'
% Sys = The control system for observer design
% ContrSys = Controller parameters (ContrSys.G1,ContrSys.G2,ContrSys.K)
%
% gammafun is defined for a scalar argument 'w' to obtain values at 1i*w

A = Sys.A;
B = Sys.B;
C = Sys.C;
% Cm = Sys.Cm;
D = Sys.D;
% Dm = Sys.Dm;

dimX = size(A,1);
dimY = size(C,1);
dimU = size(B,2);
% dimYm = size(Cm,1);

% if ~isempty(Cm)
%   L1m = L1(:,(dimY+1):end);
%   L1 = L1(:,1:dimY);
% end

if dimY ~= 1 || dimU ~= 1
  error('The system is not SISO, controller design cannot be completed.')
end

q = length(freqs);

dimZ = IMdim(freqs,dimY);

% Step 1: Compute stabilizing operators K_S and L1 so that A+B*K_S and A+L1*C
% are exponentially stable. These may also be given explicitly.
if (nargin < 3) || (isempty(K_S))
  K_S = -lqr(full(A),B,eye(size(A)),1*eye(dimU),zeros(dimX,dimU));
end
if (nargin < 4) || (isempty(L1))
  L1 = -lqr(full(A.'),C.',eye(size(A)),0.1*eye(dimY),zeros(dimX,dimY)).';
end

PKappr = @(s) (C+D*K_S)*((s*speye(dimX)-(A+B*K_S))\B)+D;
gammak = zeros(length(freqs),1);
for ind = 1:q
  PKval = PKappr(1i*freqs(ind));
  if cond(PKval)>1e6
    warning(['The value P_K(iw_k) for k=' num2str(ind) ' is nearly zero!'])
  end
  if nargin >= 5
    gammak(ind) = gammafun(freqs(ind));
    % Check accuracy of the approximation of the gamma function
    abs(gammak(ind)-PKval\((C+D*K_S)*((1i*freqs(ind)*speye(dimX)-A-B*K_S)\L1)-eye(dimY)))
%     [real(gammak(ind)-PKval\((C+D*K_S)*((1i*freqs(ind)*speye(dimX)-A-B*K_S)\L1)-eye(dimY))) imag(gammak(ind)-PKval\((C+D*K_S)*((1i*freqs(ind)*speye(dimX)-A-B*K_S)\L1)-eye(dimY)))]
%     gammak(ind)
  else
    gammak(ind) = PKval\((C+D*K_S)*((1i*freqs(ind)*speye(dimX)-A-B*K_S)\L1)-eye(dimY));
  end
end


% Construct the internal model
G1 = ConstrIM(freqs,dimY);

if freqs(1)==0
  F0 = [eye(dimY),repmat([eye(dimY),zeros(dimY)],1,q-1)];
  nzfreqs = freqs(2:end);
  Lambda = [gammak(1), zeros(dimY,2*(q-1))];
  gammaknz = gammak(2:end);
else
  F0 = repmat([eye(dimY),zeros(dimY)],1,q);
  nzfreqs = freqs;
  Lambda = zeros(dimY,2*q);
  gammaknz = gammak;
end

Lambdanz = zeros(dimY,2*length(nzfreqs));
for ind = 1:length(nzfreqs)
  
  indran = 2*(ind-1)+(1:2);
  Lambdanz(:,indran) = [real(gammaknz(ind)) imag(gammaknz(ind))];
  
end


Lambda(:,(end-(2*length(nzfreqs))+1):end) = Lambdanz;


% if freqs(1)==0
%   L2 = -place(G1',F0',[linspace(-5,-4,length(freqs)),linspace(-3.9,-3,length(nzfreqs))])';
% else
% L2 = -place(G1',F0',[linspace(-5,-4,length(freqs)),linspace(-3.9,-3,length(freqs))])';
% end
L2 = -lqr(G1.',F0.',eye(dimZ),4*eye(dimY)).';

% eig(G1+L2*F0)



% Alternative, choose [K1 K2] to stabilize ([A 0;G2*C G1],[B;G2*D])
% Kfull = -lqr([G1 G2*C;zeros(dimX,dimZ) A],[G2*D;B],blkdiag(200*eye(dimZ),.01*eye(dimX)),0.01*dimU,zeros(dimZ+dimX,dimU));
% K1 = Kfull(:,1:dimZ);
% K2 = Kfull(:,(dimZ+1):end);

% max(real(eig([G1 G2*C;zeros(dimX,dimZ) A]+[G2*D;B]*Kfull)))

ContrSys.G1 = [A+B*K_S+L1*(C+D*K_S), (B+L1*D)*Lambda; L2*(C+D*K_S), G1+L2*(F0+D*Lambda)];
ContrSys.G2 = [-L1;-L2];
ContrSys.K = [K_S, Lambda];


% L1f = [L1,L1m];
% L2f = [L2,zeros(dimZ,dimYm)];
% 
% ContrSys.G1 = [A+B*K_S+L1f*([C;Cm]+[D;Dm]*K_S), (B+L1f*[D;Dm])*Lambda; L2*(C+D*K_S), G1+L2*(F0+D*Lambda)];
% ContrSys.G2 = [-L1f;-L2f];
% ContrSys.K = [K_S, Lambda];


if issparse(A)
  ContrSys.G1 = sparse(ContrSys.G1);
end



% % Temp
% % Check validity of the controller!
% 
% Gf1 = ContrSys.G1;
% Gf2 = ContrSys.G2;
% Kf = ContrSys.K;
% Ae = [A B*Kf;Gf2*C Gf1+Gf2*D*Kf];
% Qe = [-eye(dimX) zeros(dimX,size(Gf1,2));H zeros(size(G1,1),size(Gf1,2)); -eye(dimX) zeros(dimX,size(G1,1)) eye(dimX)];


% for ind = 1:length(nzfreqs)
%   indran = zoffset+(ind-1)*2*dimY+(1:(2*dimY));
% 
%   ContrSys.G1(indran,indran) = nzfreqs(ind)*[zeros(dimY) eye(dimY);-eye(dimY) zeros(dimY)];
%   
% %   Ppi = pinv(Pvals{ind});
%   Ppi = negsqrt(Pvals{ind}); % experimental: "optimal" choice of K_0?
%   ContrSys.K(:,indran) = [real(Ppi) imag(Ppi)];
% end
% 
% if freqs(1)==0
%   ContrSys.G2 = [-eye(dimY);repmat([-eye(dimY);zeros(dimY)],length(nzfreqs),1)];
% else
%   ContrSys.G2 = repmat([-eye(dimY);zeros(dimY)],q,1);
% end
% 
% 
% 
% if length(epsgain) == 1
%   
%   ContrSys.K = epsgain*ContrSys.K;
%   return
%   
% elseif length(epsgain) == 2
%   ee_cand = linspace(epsgain(1),epsgain(2),20);
% else
%   ee_cand = epsgain;
% end
% 
% % A very crude optimization for the parameter eps>0!
% stab_margin_old = 0;
% marg_tol = 5e-4;
% 
% allmargs = zeros(size(ee_cand));
% 
% K0 = ContrSys.K;
% 
% for ind = 1:length(ee_cand)
%   ee = ee_cand(ind);
%   K = ee*K0;
%   
%   Ae = [Sys.A Sys.B*K;ContrSys.G2*Sys.C ContrSys.G1+ContrSys.G2*Sys.D*K];
%   
%   stab_margin = abs(max(real(eig(full(Ae)))));
%   
%   allmargs(ind) = stab_margin;
%   if stab_margin<stab_margin_old+marg_tol
%     
%     epsgain = ee_cand(ind-1);
%     ContrSys.K = epsgain*K0;
%     
%     break
%   else
%     stab_margin_old = stab_margin;
%     epsgain = ee;
%   end
% 
% end
% 
% function Amsq = negsqrt(A)
% % Compute "A^{-1/2}" for a matrix A
% 
% [U,S,V] = svd(A);
% Amsq = V*diag(1./sqrt(diag(S)))*U';
% 
