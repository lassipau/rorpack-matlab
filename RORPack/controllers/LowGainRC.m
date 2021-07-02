function [ContrSys,epsgain] = LowGainRC(freqsReal,Pvals,epsgain,Sys)
% Construct a low-gain simple controller for stable systems, in real form
% Inputs:
%   freqsReal : [1xN double] Frequencies to be included in the controller,
%   only real nonnegative frequencies, if zero frequency is included,
%   it's the first element in the vector
%
%   Pvals : [cell array] Values (or approximations of them) of the values
%   of the transfer function of the system on the frequencies 'freqs'
%
%   epsgain : [1x1 double/1x2 double] The value of the low-gain parameter
%   $\eps>0$, can alternatively be a vector of length 2 providing
%   minimal and maximal values for $\eps$
%
%   Sys : [struct with fields A,B,C,D] The controlled system
%
% Outputs:
%   ContrSys : [struct with fields G1,G2,K] Low-Gain Robust Controller
%   for the stable linear system Sys
%
%   epsgain : [double] The (roughly) optimal value of $\eps$

if max(real(eig(full(Sys.A))))>=0
  error('The system is unstable, the low-gain controller design cannot be completed.')
end


dimY = size(Pvals{1},1);
dimU = size(Pvals{1},2);

[G1,G2tmp] = ConstrIM(freqsReal,dimY);

ContrSys.G1 = G1;
ContrSys.G2 = -G2tmp;
dimZ = size(ContrSys.G1,1);
  
ContrSys.K = zeros(dimU,dimZ);

if freqsReal(1)==0
  zoffset = dimY; 
  ContrSys.K(:,1:dimY) = pinv(Pvals{1});
  nzfreqs = freqsReal(2:end);
else
  zoffset = 0;
  nzfreqs = freqsReal;
end

for ind = 1:length(nzfreqs)
  indran = zoffset+(ind-1)*2*dimY+(1:(2*dimY));
  Ppi = pinv(Pvals{ind});
  ContrSys.K(:,indran) = [real(Ppi), imag(Ppi)];
end


if length(epsgain) == 1
  
  ContrSys.K = epsgain*ContrSys.K;
  return
  
elseif length(epsgain) == 2
  ee_cand = linspace(epsgain(1),epsgain(2),20);
else
  ee_cand = epsgain;
end

% A very crude optimization for the parameter eps>0!
stab_margin_old = -1;
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