function [ContrSys,epsgain] = ConstrContrLG(freqs,Pvals,epsgain,Sys)
% ContrSys = ConstrContrLG(freqs,dimY,Pvals)
%
% Construct a low-gain simple controller for stable systems
% freqs = Frequencies to be included in the controller
% Pvals = [cell array] Values (or approximations of them) of the values of the transfer
% function of the system on the frequencies 'freqs'
% ContrSys = Controller parameters (ContrSys.G1,ContrSys.G2,ContrSys.K)

if max(real(eig(full(Sys.A))))>=0
  error('The system is unstable, the low-gain controller design cannot be completed.')
end


dimY = size(Pvals{1},1);
dimU = size(Pvals{1},2);
q = length(freqs);

ContrSys.G1 = zeros(q*dimY);
ContrSys.K = zeros(dimU,size(ContrSys.G1,1));

for ind = 1:q
  indran = (ind-1)*dimY+(1:dimY);

  ContrSys.G1(indran,indran) = freqs(ind)*eye(dimY);
  
  ContrSys.K(:,indran) = pinv(Pvals{ind});
  % ContrSys.K(:,indran) = negsqrt(Pvals{ind}); % experimental: "optimal" choice of K_0?
end


ContrSys.G2 = repmat(-eye(dimY),q,1);

if length(epsgain) == 1
  
  ContrSys.K = epsgain*ContrSys.K;
  return
  
elseif length(epsgain) == 2
  ee_cand = linspace(epsgain(1),epsgain(2),20);
else
  ee_cand = epsgain;
end

% A very crude optimization for the parameter eps>0!
stab_margin_old = 0;
marg_tol = 5e-3;

allmargs = zeros(size(ee_cand));

K0 = ContrSys.K;

for ind = 1:length(ee_cand)
ee = ee_cand(ind);
K = ee*K0;

Ae = [Sys.A Sys.B*K;ContrSys.G2*Sys.C ContrSys.G1+ContrSys.G2*Sys.D*K];

stab_margin = abs(max(real(eig(full(Ae)))));

allmargs(ind) = stab_margin;
if stab_margin<stab_margin_old+marg_tol
    
    epsgain = ee_cand(ind-1);
    ContrSys.K = epsgain*K0;

    break
else
    stab_margin_old = stab_margin;
    epsgain = ee;
end

end

function Amsq = negsqrt(A)
% Compute "A^{-1/2}" for a matrix A

[U,S,V] = svd(A);
Amsq = V*diag(1./sqrt(diag(S)))*U';
