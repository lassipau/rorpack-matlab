function K = lqr_RP(A,B,C,R)
% function K = lqr_RP(A,B,C,R)
% 
% Use LQR design to compute K such that A-B*K is Hurwitz. 
% The matrix K is given by K=R\B'X, where X is the positive solution of the 
% algebraic Riccati equation
% 
% X*A + A.'*X - X*B*(R\B.')*X + C.'*C = 0
% 
% The solution uses 'icare' for Matlab versions R2019a and newer, and
% 'care' for earlier Matlab versions.
%
% Copyright (C) 2020 by Lassi Paunonen (lassi.paunonen@tuni.fi)
% Licensed under GNU GPLv3 (see LICENSE.txt).


dimX = size(A,1);
dimU = size(B,2);
dimY = size(C,1);

if verLessThan('matlab','R2019a')
    
    warning('Matlab version prior to R2019a: Using "care" instead of "icare" (may be less accurate or slower).')
    
    B_ext = [B,zeros(dimX,dimY)];
    Ss_ext = [zeros(dimX,dimU),C'];
    R_ext = blkdiag(R,-eye(dimY));
    [~,evals,K_ext] = care(A,B_ext,0,R_ext,Ss_ext);
    K = K_ext(1:dimU,:);
    
    if max(real(evals))>= 0
        error('Stabilization of the pair (A,B) failed!')
    end
    
else
   
    B_ext = [B,zeros(dimX,dimY)];
    Ss_ext = [zeros(dimX,dimU),C'];
    R_ext = blkdiag(R,-eye(dimY));
    [~,K_ext,evals] = icare(A,B_ext,0,R_ext,Ss_ext,[],[]);
    K = K_ext(1:dimU,:);
    
    if max(real(evals))>= 0
        error('Stabilization of the pair (A,B) failed!')
    end
    
    
end





