function ContrSys = StabilizingObserverBased(Sys,K,L)
% function ContrSys = StabilizingObserverBased(Sys,K,L)
%
% Construct an stabilizing error feedback controller based on Luenberger 
% observer-based design.
%
% Inputs:
%   K : [M1xN1 matrix] A stabilizing feedback gain for the pair (A,B) (such 
%                      that A+BK is stable).
%
%   L : [M2xN2 matrix] A stabilizing output injection gain for the pair
%                      (C,A) (such that A+LC is stable).
%
% Outputs:
%   ContrSys : [struct with fields G1,G2,K] observer-based output feedback 
%   controller stabilizing the closed-loop system.


ContrSys.G1 = Sys.A+L*Sys.C+Sys.B*K;
ContrSys.G2 = -L;
ContrSys.K = K;

% TODO: Add a "controller disturbance input" ContrSys.G2d so that the
% "wdist(t)" signal can be used as an external control input to the
% stabilized system. This requires defining Bd=B in the model itself.


if issparse(Sys.A)
  ContrSys.G1 = sparse(ContrSys.G1);
end
