function [ContrSys,K21] = ObserverBasedRC(freqsReal,Sys,K21,L,IMstabtype,IMstabmarg)
% function [ContrSys,K21] = ObserverBasedRC(freqsReal,Sys,K21,L,IMstabtype,IMstabmarg)
%
% Construct an Observer-Based Robust Controller for systems with
% the same number of inputs and outputs. The construction follows the 
% article "Controller Design for Robust Output Regulation of Regular Linear 
% Systems" by L. Paunonen, IEEE TAC, 2016 (Section VI) where it was
% presented for regular linear systems (generalising earlier results by
% T. Hämäläinen and S. Pohjolainen, SICON, 2010, and E. Immonen, SICON, 2007.
%
% Inputs:
%   freqsReal : [1xN double] Frequencies to be included in the controller,
%   only real nonnegative frequencies, if zero frequency is included,
%   it's the first element in the vector
%
%   Sys : [struct with fields A,B,C,D] The controlled system
%
%   K21 : [M1xN1 matrix] The matrix K21 which should be chosen
%   so that A + B * K21 is stable. If 'IMstabtype' is 'full_K', then this
%   parameter is instead an operator K chosen so that A_s + B_s * K is
%   exponentially stable, where A_s = [G1,G2*C;0,G1] and B_s=[G2*D;B].
%
%   L : [M2xN2 matrix] The matrix L which should be chosen
%   so that A + L * C is stable.
%
%   IMstabtype : [string] Stabilization of the internal model
%   using either 'LQR', 'poleplacement', or 'full_K'.
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

q = length(freqsReal);


if ~isequal(IMstabtype,'full_K')
    PKappr = @(s) (C+D*K21)*((s*eye(dimX)-(A+B*K21))\B)+D;
    for ind = 1:q
        if cond(PKappr(1i*freqsReal(ind)))>1e6
            warning(['The matrix P_K(iw_k) for k=' num2str(ind) ' is nearly singular!'])
        end
    end
end

% Construct the internal model
[G1,G2] = ConstrIM(freqsReal,dimY);

dimZ = size(G1,1);

if isequal(IMstabtype,'full_K')
    % The full operator K is given in the place of parameter "K21"
    
    % Check that the dimensions are correct
    if ~isequal(size(K21),[dimU,dimZ+dimX])
        error('The dimensions of the input parameter K21 are not correct!')
    end
    
    K1 = K21(:,1:dimZ);
    K2 = K21(:,(dimZ+1):end);
    
else
    
    if dimY ~= dimU
      error('The system has an unequal number of inputs and outputs, the observer-based controller design cannot be completed (in this form).')
    end

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
            target_eigs = target_eigs(:);
            K1 = -place(G1,B1,target_eigs);
        else
            target_eigs = ones(2*length(freqsReal),1)*linspace(-1.1*IMstabmarg,-IMstabmarg,dimY)...
                +1i*reshape([freqsReal(end:-1:1),-freqsReal],2*length(freqsReal),1)*ones(1,dimY);
            target_eigs = target_eigs(:);
            K1 = -place(G1,B1,target_eigs);
        end
    else
        error('Unknown stabilization type for the Internal Model')
    end
    
    K2 = K21+K1*H;
end


% IN PROGRESS: An alternative option for stabilization: If K21 not given, 
% could use LQR to construct [K1 K2] as a stabilizing feedback gain for the 
% pair ([G1,G2*C;0,G1],[G2*D;B])
% Kfull = -lqr([G1 G2*C;zeros(dimX,dimZ) A],[G2*D;B],blkdiag(200*eye(dimZ),.01*eye(dimX)),0.01*eye(dimU),zeros(dimZ+dimX,dimU));
% K1 = Kfull(:,1:dimZ);
% K2 = Kfull(:,(dimZ+1):end);
% max(real(eig([G1 G2*C;zeros(dimX,dimZ) A]+[G2*D;B]*Kfull)))


ContrSys.G1 = [G1 zeros(dimZ,dimX); (B+L*D)*K1 A+B*K2+L*(C+D*K2)];
ContrSys.G2 = [G2;-L];
ContrSys.K = [K1, K2];
ContrSys.Dc = zeros(dimU,dimY);

if issparse(A)
  ContrSys.G1 = sparse(ContrSys.G1);
end



