function [ContrSys,epsgain] = PassiveRC(freqs,Pvals,epsgain,Sys)
% ContrSys = PassiveRC(freqs,dimY,Pvals)
%
% Construct a passive simple controller for stable systems, in real form
% freqs = Frequencies to be included in the controller, only real nonnegative
% frequencies, if zero frequency is included, it's the first element in the
% vector
% Pvals = [cell array] Values (or approximations of them) of the values of the transfer
% function of the system on the frequencies 'freqs'
% ContrSys = Controller parameters (ContrSys.G1,ContrSys.G2,ContrSys.K)

%if max(real(eig(full(Sys.A))))>=0
%  error('The system is unstable, the low-gain controller design cannot be completed.')
%end


dimY = size(Pvals{1},1);
dimU = size(Pvals{1},2);
q = length(freqs);

dimZ = IMdim(freqs,dimY);
  
  
ContrSys.G1 = zeros(dimZ);
ContrSys.K = zeros(dimU,dimZ);

% To add: weights for the components in B_c, based on values of the
% transfer function at the frequencies


if freqs(1)==0
  zoffset = dimY; 
  ContrSys.K(:,1:dimY) = eye(dimY);
%   ContrSys.K(:,1:dimY) = negsqrt(Pvals{1}); % experimental: "optimal" choice of K_0?
  nzfreqs = freqs(2:end);
else
  zoffset = 0;
  nzfreqs = freqs;
end

for ind = 1:length(nzfreqs)
  indran = zoffset+(ind-1)*2*dimY+(1:(2*dimY));

  ContrSys.G1(indran,indran) = nzfreqs(ind)*[zeros(dimY) eye(dimY);-eye(dimY) zeros(dimY)];
  
%   Ppi = pinv(Pvals{ind});
%   Ppi = negsqrt(Pvals{ind}); % experimental: "optimal" choice of K_0?
%   ContrSys.K(:,indran) = [real(Ppi) imag(Ppi)];
end

if freqs(1)==0
  ContrSys.G2 = [-eye(dimY);repmat([-eye(dimY);zeros(dimY)],length(nzfreqs),1)];
else
  ContrSys.G2 = repmat([-eye(dimY);zeros(dimY)],q,1);
end

ContrSys.K = -(ContrSys.G2).';

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
marg_tol = 5e-4;

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

