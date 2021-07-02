function ContrSys = ObserverBasedROMRC(freqsReal,SysApprox,alpha1,alpha2,R1,R2,Q0,Q1,Q2,ROMorder)
% Construct a reduced order observer-based robust controller for systems
% with the same number of inputs and outputs.
% The frequencies are assumed to be conjugate pairs, and the internal
% model is in real form
% Inputs:
%   freqsReal : [1xN double] Frequencies to be included in the controller,
%   only real nonnegative frequencies, if zero frequency is included,
%   it's the first element in the vector.
%
%   SysApprox : The Galerkin approximation (A^N,B^N,C^N,D) of the control
%   system for the controller design. This is a struct with fields
%   SysApprox.AN, SysApprox.BN, SysApprox.CN, and SysApprox.D
%
%   IMstabmarg : [double] intended stability margin of the closed-loop system
%
%   ROMorder : order of the reduced-order observer in the controller. The
%            model reduction step can be skipped by setting 'ROMorder=NaN'
%
% In this version, the parameters in the LQR/LQG design Q are chosen to be
% identity matrices.
%
% Output:
%   ContrSys : [struct with fields G1,G2,K] reduced order
%   observer-based robust controller
%
% Copyright (C) 2020 by Lassi Paunonen (lassi.paunonen@tuni.fi)
% Licensed under GNU GPLv3 (see LICENSE.txt).

AN = SysApprox.AN;
BN = SysApprox.BN;
CN = SysApprox.CN;
D = SysApprox.D;

dimX = size(AN,1);
dimY = size(CN,1);
dimU = size(BN,2);

% Check the consistency of the controller design parameters

alpha1_check = isreal(alpha1) && isscalar(alpha1) && alpha1>0;
alpha2_check = isreal(alpha2) && isscalar(alpha2) && alpha2>0;

R1_check = ishermitian(R1) && all(eig(R1)>0);
R2_check = ishermitian(R2) && all(eig(R2)>0);

if ~alpha1_check || ~alpha2_check
    error('"alpha1" and "alpha2" need to be positive!')
elseif ~R1_check || ~R2_check
    error('"R1" and "R2" need to be positive definite')
end

% Construct the internal model
[G1,G2] = ConstrIM(freqsReal,dimY);
dimZ0 = size(G1,1);


% If R1, R2, Q0, Q1, or Q2 have scalar values, these are interpreted as
% "scalar times identity".
if dimY>1 && isscalar(R1)
    R1 = R1*eye(dimY);
end
if dimY>1 && isscalar(R2)
    R2 = R2*eye(dimY);
end
if dimZ0>1 && isscalar(Q0)
    Q0 = Q0*eye(dimZ0);
end
if dimX>1 && isscalar(Q1)
    Q1 = Q1*eye(dimX);
end
if dimX>1 && isscalar(Q2)
    Q2 = Q2*eye(dimX);
end


% Check the consistency of the dimensions of R1, R2, Q0, Q1, and Q2
if ~isequal(size(R1),[dimY,dimY])
    error('Dimensions of "R1" are incorrect! (should be [dimY,dimY]).')
end
if ~isequal(size(R2),[dimU,dimU])
    error('Dimensions of "R1" are incorrect! (should be [dimU,dimU]).')
end
if ~isequal(size(Q0,2),dimZ0)
    error('Dimensions of "Q0" are incorrect!')
end
if ~isequal(size(Q1,1),dimX)
    error('Dimensions of "Q1" are incorrect!')
end
if ~isequal(size(Q2,2),dimX)
    error('Dimensions of "Q2" are incorrect!')
end


% Form the extended system (As,Bs)
As = [G1,G2*CN;zeros(dimX,dimZ0),AN];
Bs = [G2*D;BN];

Qs = blkdiag(Q0,Q2);


% Stabilize the pairs (CN,AN+alpha1) and (As+alpha2,Bs) using LQR/LQG design
if verLessThan('matlab','R2019a')
    
    warning('Matlab version prior to R2019a: Using "care" instead of "icare" (may be less accurate or slower).')
    
    % Stabilize the pair (CN,AN+alpha1)
    tic
    B_ext = [CN',zeros(dimX)];
    S_ext = [zeros(dimX,dimY),Q1];
    R1_ext = blkdiag(R1,-eye(dimX));
    [~,evals,L_ext_adj] = care((AN+alpha1*eye(dimX))',B_ext,0,R1_ext,S_ext);
    fprintf(['Stabilization of the pair (CN,AN+alpha1) took ' num2str(round(toc,1)) ' seconds.\n'])
    L = -L_ext_adj(1:dimY,:)';
    
    if max(real(evals))>= 0
        error('Stabilization of the pair (CN,AN+alpha1) failed!')
    end
    
    % Stabilize the pair (As+alpha2,Bs)
    tic
    Bs_ext = [Bs,zeros(dimZ0+dimX)];
    Ss_ext = [zeros(dimZ0+dimX,dimY),Qs'];
    R2_ext = blkdiag(R2,-eye(dimZ0+dimX));
    [~,evals,K_ext] = care(As+alpha2*eye(dimZ0+dimX),Bs_ext,0,R2_ext,Ss_ext);
    K = -K_ext(1:dimU,:);
    fprintf(['Stabilization of the pair (As+alpha2,Bs) took ' num2str(round(toc,1)) ' seconds.\n'])
    
    if max(real(evals))>= 0
        error('Stabilization of the pair (As+alpha2,Bs) failed!')
    end
    
else
    
    % Stabilize the pair (CN,AN+alpha1)
    tic
    B_ext = [CN',zeros(dimX)];
    S_ext = [zeros(dimX,dimY),Q1];
    R1_ext = blkdiag(R1,-eye(dimX));
    [~,L_ext_adj,evals] = icare((AN+alpha1*eye(dimX))',B_ext,0,R1_ext,S_ext,[],[]);
    L = -L_ext_adj(1:dimY,:)';
    fprintf(['Stabilization of the pair (CN,AN+alpha1) took ' num2str(round(toc,1)) ' seconds.\n'])
    
    if max(real(evals))>= 0
        error('Stabilization of the closed-loop system failed!')
    end
    
    % Stabilize the pair (As+alpha2,Bs)
    tic
    Bs_ext = [Bs,zeros(dimZ0+dimX)];
    Ss_ext = [zeros(dimZ0+dimX,dimY),Qs'];
    R2_ext = blkdiag(R2,-eye(dimZ0+dimX));
    [~,K_ext,evals] = icare(As+alpha2*eye(dimZ0+dimX),Bs_ext,0,R2_ext,Ss_ext,[],[]);
    K = -K_ext(1:dimU,:);
    fprintf(['Stabilization of the pair (As+alpha2,Bs) took ' num2str(round(toc,1)) ' seconds.\n'])
    
    if max(real(evals))>= 0
        error('Stabilization of the pair (As+alpha2,Bs) failed!')
    end
    
    
end

% Decompose the control gain K into K=[K1N,K2N]
K1N = K(:,1:dimZ0);
K2N = K(:,(dimZ0+1):end);


% Complete the model reduction step of the controller design. If the model
% reduction fails (for example due to too high reduction order), the
% controller design is by default completed without the model reduction.
if ~isnan(ROMorder)
    try
        tic
        rsys = balred(ss(AN+L*CN,[BN+L*D,L],K2N,zeros(dimU,dimU+dimY)),ROMorder);
        fprintf(['Model reduction step took ' num2str(round(toc,1)) ' seconds.\n'])
        
        ALr = rsys.A;
        Br_full = rsys.B;
        BLr = Br_full(:,1:dimU);
        Lr = Br_full(:,(dimU+1):end);
        K2r = rsys.C;
    catch
        warning(['Model reduction step failed! Modify "ROMorder", or check that "balred" is available. ', ...
            'Proceeding without model reduction (with option "ROMorder=NaN")']);
    end
end

if isnan(ROMorder)
    fprintf('Constructing the controller without model reduction.\n')
    
    ALr = AN+L*CN;
    BLr = BN+L*D;
    Lr = L;
    K2r = K2N;
end


% Construct the controller matrices (G1,G2,K).
ContrSys.G1 = [G1,zeros(dimZ0,size(ALr,1));BLr*K1N,ALr+BLr*K2r];
ContrSys.G2 = [G2;-Lr];
ContrSys.K = [K1N, K2r];








