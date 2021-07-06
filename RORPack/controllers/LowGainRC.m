function [ContrSys,epsgain] = LowGainRC(freqsReal,Pvals,epsgain,Sys,Dc)
% function [ContrSys,epsgain] = LowGainRC(freqsReal,Pvals,epsgain,Sys)
%
% Construct a low-gain simple controller for stable systems. The parameters
% are chosen as in the article "Controller Design for Robust Output
% Regulation of Regular Linear Systems" by L. Paunonen IEEE TAC 2016
% (Section IV).
%
% Inputs:
%   freqsReal : [1xN double] Frequencies to be included in the controller,
%   only real nonnegative frequencies, if zero frequency is included,
%   it's the first element in the vector
%
%   Pvals : [cell array] Values (or approximations of them) of the values
%   of the transfer function (A,B,C,D) UNDER OUTPUT FEEDBACK u(t)=Dc*y(t) 
%   at the complex frequencies in 'freqsReal', i.e., P_{Dc}(i*w_k). If the
%   transfer function P(s) of the system (A,B,C,D) is well-defined at
%   i*w_k, then these values can be computed with the formula 
%   P_{Dc}(i*w_k) = (I-P(i*w_k)*Dc)\P(i*w_k). If Dc=0, then in particular
%   P_{Dc}(i*w_k) = P(i*w_k).
%
%   epsgain : [1x1 double/1x2 double] The value of the gain parameter
%   $\eps>0$, which scales both the input operator G2 and the output 
%   operator K. The parameter can alternatively be a vector of length 2 
%   providing minimal and maximal values for $\eps$ between which this
%   parameter is optimized (crudely).
%
%   Sys : [struct with fields A,B,C,D] The controlled system (used only for
%   optimizing the gain parameter if 'epsgain' is a vector with two
%   elements. If 'epsgain' is a scalar, can be omitted or replaced with '[]'.
%
%   Dc : [NxN double/1x1 double] (optional) An optional feedthrough term 
%   for the controller. The closed-loop structure implies that a nonzero
%   feedthrough term with 'Dc' is equivalent to a 'prestabilizing' output
%   feedback to the system (plus an additional disturbance term). This term 
%   is required if the system is initially unstable, and the output 
%   feedback u(t)=Dc*y(t) is required to lead to an exponentially stable 
%   system. 
%
% Outputs:
%   ContrSys : [struct with fields G1,G2,K] Low-Gain Robust Controller
%   for the stable linear system Sys
%
%   epsgain : [double] The (roughly) optimal value of $\eps$


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

% Define the feedthrough term of the controller. If Dc is not given as a
% parameter, we define Dc=0.
if nargin > 4
    if ~isequal(size(Dc),[dimU,dimY]) && ~isequal(Dc,0)
        error('The controller feedthrough Dc has wrong dimensions (should be dimU x dimY)!')
    elseif isequal(Dc,0)
        ContrSys.Dc = zeros(dimU,dimY);
    else
        ContrSys.Dc = Dc;
    end
    
    % Check that the controller feedthrough term stabilizes the original
    % system.
    SysFB = SysOutputFeedback(Sys,ContrSys.Dc);
    if max(real(eig(full(SysFB.A))))>=0
        error('The output feedback u(t)=Dc*y(t) does not stabilize the system!')
    end
else
    ContrSys.Dc = zeros(dimU,dimY);
    % Check that the system is exponentially stable.
    if max(real(eig(full(Sys.A))))>=0
        error('The system is unstable, the passive controller design cannot be completed.')
    end
end


% Determine the gain parameter either by using the given 'epsgain' (if the
% parameter is a scalar), or by optimizing over the given range of values.
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
marg_tol = 5e-5;

allmargs = zeros(size(ee_cand));

K0 = ContrSys.K;

for ind = 1:length(ee_cand)
  ee = ee_cand(ind);
  ContrSys.K = ee*K0;
  
  CLSys_tmp = ConstrCLSys(Sys,ContrSys);
  
  stab_margin = abs(max(real(eig(full(CLSys_tmp.Ae)))));
  
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