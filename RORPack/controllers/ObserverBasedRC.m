function [ContrSys,K21] = ObserverBasedRC(freqsReal,Sys,K21,L,IMstabtype,IMstabmarg)
% Construct an observer-based robust controller for systems with
% the same number of inputs and outputs. 
% Inputs:
%   freqsReal : [1xN double] Frequencies to be included in the controller,
%   only real nonnegative frequencies, if zero frequency is included,
%   it's the first element in the vector
%
%   Sys : [struct with fields A,B,C,D] The controlled system
%
%   K21 : [M1xN1 matrix] The matrix K21 which should be chosen
%   so that A + B * K21 is stable.
%
%   L : [M2xN2 matrix] The matrix L which should be chosen
%   so that A + L * C is stable.
%
%   IMstabtype : [string] Stabilization of the internal model
%   using either 'LQR' or 'poleplacement'
%
%   IMstabmarg : [double] The desired stability margin
%   for the internal model.
%
% Outputs:
%   ContrSys : [struct with fields G1,G2,K] observer-based robust controller
%   for the stable linear system Sys

A = Sys.A;
B = Sys.B;
C = Sys.C;
D = Sys.D;

dimX = size(A,1);
dimY = size(C,1);
dimU = size(B,2);

if dimY ~= dimU
  error('The system has an unequal number of inputs and outputs, the observer-based controller design cannot be completed.')
end

q = length(freqsReal);


%if freqs(1)==0, dimZ = dimY*(2*q-1); else dimZ = dimY*2*q; end

%B1 = zeros(dimZ,dimY);
%if freqs(1)==0
%  offset = 1;
%  B1(1:dimY,:) = PKvals;
%end


PKappr = @(s) (C+D*K21)*((s*eye(dimX)-(A+B*K21))\B)+D;
for ind = 1:q
  if cond(PKappr(1i*freqsReal(ind)))>1e6
    warning(['The matrix P_K(iw_k) for k=' num2str(ind) ' is nearly singular!'])
  end
end


% Construct the internal model
[G1,G2] = ConstrIM(freqsReal,dimY);

dimZ = size(G1,1);

% Find H as the solution of G1*H=H*(A+B*K21)+G2*(C+D*K21) and define B1
H = sylvester(G1,-(A+B*K21),G2*(C+D*K21));
B1 = H*B+G2*D;

% Stabilization of the internal model, choose K1 so that G1+B1*K1 is
% exponentially stable
if isequal(IMstabtype,'LQR')
  K1 = -lqr(G1+IMstabmarg*eye(dimZ),B1,100*eye(dimZ),0.001*eye(dimU),zeros(dimZ,dimU));
elseif isequal(IMstabtype,'poleplacement')
  if freqsReal(1)==0
    target_eigs = ones(2*length(freqsReal)-1,1)*linspace(-1.1*IMstabmarg,-IMstabmarg,dimY)...
                    +1i*reshape([freqsReal(end:-1:1),-freqsReal(2:end)],2*length(freqsReal)-1,1)*ones(1,dimY);
%     target_eigs
%     target_eigs = [linspace(-7,-5,length(freqs)*dimY),linspace(-4.9,-3,length(nzfreqs)*dimY)];
    target_eigs = target_eigs(:);
    K1 = -place(G1,B1,target_eigs);
  else
    target_eigs = ones(2*length(freqsReal),1)*linspace(-1.1*IMstabmarg,-IMstabmarg,dimY)...
                    +1i*reshape([freqsReal(end:-1:1),-freqsReal],2*length(freqsReal),1)*ones(1,dimY);
%     target_eigs = [linspace(-5,-4,length(freqs)*dimY),linspace(-3.9,-3,length(freqs)*dimY)];
    target_eigs = target_eigs(:);
    K1 = -place(G1,B1,target_eigs);
  end
else
  error('Unknown stabilization type for the Internal Model')
end
% Add: scaling to K1 or G2 through an additional gain parameter!

K2 = K21+K1*H;


% Alternative, choose [K1 K2] to stabilize ([A 0;G2*C G1],[B;G2*D])
% Kfull = -lqr([G1 G2*C;zeros(dimX,dimZ) A],[G2*D;B],blkdiag(200*eye(dimZ),.01*eye(dimX)),0.01*eye(dimU),zeros(dimZ+dimX,dimU));
% K1 = Kfull(:,1:dimZ);
% K2 = Kfull(:,(dimZ+1):end);

% max(real(eig([G1 G2*C;zeros(dimX,dimZ) A]+[G2*D;B]*Kfull)))


ContrSys.G1 = [G1 zeros(dimZ,dimX); (B+L*D)*K1 A+B*K2+L*(C+D*K2)];
ContrSys.G2 = [G2;-L];
ContrSys.K = [K1, K2];


if issparse(A)
  ContrSys.G1 = sparse(ContrSys.G1);
end



