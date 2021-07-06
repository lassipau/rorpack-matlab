function [ContrSys,epsgain] = PassiveRC(freqsReal,dimY,epsgain,Sys,Dc)
% function [ContrSys,epsgain] = PassiveRC(freqsReal,dimY,epsgain,Sys,Dc)
%
% Construct a Passive Robust Controller for an impedance passive linear
% system which is either exponentially or strongly stable, or stabilizable 
% (i.e., can be stabilized with negative output feedback).
%
% Inputs:
%   freqsReal : [1xN double] Frequencies to be included in the controller,
%   only real nonnegative frequencies, if zero frequency is included,
%   it's the first element in the vector
%
%   dimY : [integer] Dimension of the output space Y, i.e., number of 
%   outputs of the system (note that by impedance passivity, the system 
%   has the same number of inputs and outputs).
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
%   feedback to the system (plus an additional disturbance term). In
%   particular, controllable unstable passive systems can be prestabilized 
%   by setting a negative definite matrix 'Dc'. This term is required if
%   the system is initially unstable, and the output feedback u(t)=Dc*y(t)
%   is required to lead to an exponentially or strongly stable system.
%   If 'Dc' is a scalar, it is interpreted as "Dc*eye(dimY)".
%
% Outputs:
%   ContrSys : [struct with fields G1,G2,K] Passive Robust Controller
%   for the stable linear system Sys
%
%   epsgain : [double] The (roughly) optimal value of $\eps$



[G1,G2tmp] = ConstrIM(freqsReal,dimY);

ContrSys.G1 = G1;
ContrSys.G2 = -G2tmp;
ContrSys.K = G2tmp.';

if nargin > 4
    if isscalar(Dc)
        ContrSys.Dc = Dc*eye(dimY);
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
    ContrSys.Dc = zeros(dimY);
    % Check that the system is exponentially stable.
    if max(real(eig(full(Sys.A))))>=0
        error('The system is unstable, the passive controller design cannot be completed.')
    end
end



if length(epsgain) == 1
  
  ContrSys.K = epsgain*ContrSys.K;
  ContrSys.G2 = epsgain*ContrSys.G2;
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

% In the passive controller design, we scale both the input and output
% operators by 'epsgain'
K0 = ContrSys.K;
G20 = ContrSys.G2;

for ind = 1:length(ee_cand)
  ee = ee_cand(ind);
  
  ContrSys.K = ee*K0;
  ContrSys.G2 = ee*G20;
  
  CLSys_tmp = ConstrCLSys(Sys,ContrSys);
  
  stab_margin = abs(max(real(eig(full(CLSys_tmp.Ae)))));
  
  allmargs(ind) = stab_margin;
  if stab_margin<stab_margin_old+marg_tol
    
    epsgain = ee_cand(ind-1);
    ContrSys.K = epsgain*K0;
    ContrSys.G2 = epsgain*G20;
    
    break
  else
    stab_margin_old = stab_margin;
    epsgain = ee;
  end

end

