function [ContrSys,K2] = DualObserverBasedRC(freqsReal,Sys,K2,L1,IMstabtype,IMstabmarg)
% function [ContrSys,K2] = DualObserverBasedRC(freqsReal,Sys,K2,L1,IMstabtype,IMstabmarg)
%
% Construct the Dual Observer-Based Robust Controller for systems with
% the same number of inputs and outputs. The design was introduced in the
% article "Controller Design for Robust Output Regulation of Regular Linear 
% Systems" by L. Paunonen, IEEE TAC, 2016 (Section V).
%
% Inputs:
%   freqsReal : [1xN double] Frequencies to be included in the controller,
%   only real nonnegative frequencies, if zero frequency is included,
%   it's the first element in the vector
%
%   Sys : [struct with fields A,B,C,D] The controlled system
%
%   K2 : [M1xN1 matrix] The matrix K2 which should be chosen
%   so that A + B * K2 is stable.
%
%   L : [M2xN2 matrix] The matrix L1 which should be chosen
%   so that A + L1 * C is stable.
%
%   IMstabtype : [string] Stabilization of the internal model
%   using either 'LQR' or 'poleplacement'
%
%   IMstabmarg : [double] The desired stability margin
%   for the internal model.
%
% Outputs:
%   ContrSys : [struct with fields G1,G2,K] dual observer-based robust controller
%   for the stable linear system Sys

A = Sys.A;
B = Sys.B;
C = Sys.C;
D = Sys.D;

dimX = size(A,1);
dimY = size(C,1);
dimU = size(B,2);

q = length(freqsReal);


PLappr = @(s) C*((s*eye(dimX)-(A+L1*C))\(B+L1*D))+D;
for ind = 1:q
  if rank(PLappr(1i*freqsReal(ind)),1e-6)<dimY
    warning(['The matrix P_L(iw_k) for l=' num2str(ind) ' is (at least nearly) non-surjective!'])
  end
  if max(max(abs(conj(PLappr(1i*freqsReal(ind)))-PLappr(-1i*freqsReal(ind)))))>1e-4
    warning('The controller design is only guaranteed to work for systems satisfying P(-iw_k)=conj(P(-iw_k))!')
  end

end


% Construct the internal model
[G1,G2] = ConstrIM(freqsReal,dimY);
dimZ = size(G1,1);


if dimY == dimU
    %   error('The system has an unequal number of inputs and outputs, the observer-based controller design cannot be completed.')
    % If the number of inputs and outputs are the same, choose K using
    % blocks of identity and zero matrices
    K1 = G2.';
else 
    % Otherwise choose parts K1_k of K1 in such a way that P_L(iw_k)K1_k are
    % invertible. These choices are based on the approximate values of the
    % transfer function P_L(s) at iw_k. Note that using approximations is valid
    % in this part of the controller construction, since the purpose is
    % only to guarantee the stabilizability of the pair (C_1,G_1).
    
    K1 = zeros(dimU,dimZ);
    
    if freqsReal(1)==0
        K1(:,1:dimY) = pinv(PLappr(0));
        ind_offset = 1;
    else
        ind_offset = 0;
    end
    
    for ind = (ind_offset+1):q
        PLpinv = pinv(PLappr(1i*freqsReal(ind)));
        K1(:,(2*(ind-1-ind_offset)+ind_offset)*dimY+(1:(2*dimY))) = [real(PLpinv),imag(PLpinv)];
    end
end



% Find H as the solution of H*G1=(A+L1*C)*H+(B+L1*D)*K and define C1
H = sylvester(-(A+L1*C),G1,(B+L1*D)*K1);
C1 = C*H+D*K1;

% Stabilization of the internal model, choose G2 so that G1+G2*C1 is
% exponentially stable
if isequal(IMstabtype,'LQR')
  G2 = conj(-lqr(conj(G1).'+IMstabmarg*eye(dimZ),conj(C1).',100*eye(dimZ),0.001*eye(dimY),zeros(dimZ,dimY))).';
elseif isequal(IMstabtype,'poleplacement')
  if freqsReal(1)==0
    target_eigs = ones(2*length(freqsReal)-1,1)*linspace(-1.1*IMstabmarg,-IMstabmarg,dimY)...
                    +1i*reshape([freqsReal(end:-1:1),-freqsReal(2:end)],2*length(freqsReal)-1,1)*ones(1,dimY);
    target_eigs = target_eigs(:);
    G2 = conj(-place(conj(G1).',conj(C1).',target_eigs)).';
  else
    target_eigs = ones(2*length(freqsReal),1)*linspace(-1.1*IMstabmarg,-IMstabmarg,dimY)...
                    +1i*reshape([freqsReal(end:-1:1),-freqsReal],2*length(freqsReal),1)*ones(1,dimY);
    target_eigs = target_eigs(:);
    G2 = conj(-place(conj(G1).',conj(C1).',target_eigs)).';
  end
else
  error('Unknown stabilization type for the Internal Model')
end

L = L1+H*G2;

% Construct the controller parameters
ContrSys.G1 = [G1 G2*(C+D*K2); zeros(dimX,dimZ) A+B*K2+L*(C+D*K2)];
ContrSys.G2 = [G2;L];
ContrSys.K = [K1, -K2];
ContrSys.Dc = zeros(dimU,dimY);

if issparse(A)
  ContrSys.G1 = sparse(ContrSys.G1);
end

