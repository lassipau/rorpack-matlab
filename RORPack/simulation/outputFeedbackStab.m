function Sys = outputFeedbackStab(Sys, K)
% Computes the linear system 
% (A+BK(I+DK)^(-1)C,B(I+KD)^(-1),(I+DK)^(-1)C,(I+DK)^(-1)D) 
% obtained from (A,B,C,D) with output feedback u(t)=Ky(t)

% K is the output feedback operator.
% Sys is the linear system.

dim_U = size(Sys.D,2);
dim_Y = size(Sys.D,1);

Q = inv(eye(dim_Y) - Sys.D*K);
Q2 = eye(dim_U) + K*(Q*Sys.D);
Sys.B = Sys.B*Q2;
Sys.C = Q*Sys.C;
Sys.D = Q*Sys.D;
Sys.A = Sys.A + Sys.B*(K*Sys.C);

% temp = sys.D*K;
% res = inv(eye(size(temp,1), size(temp,2)) - temp);
% D0 = res*sys.D;
% Dd0 = res*sys.Dd;
% Bd0 = sys.Bd + sys.B*(K*Dd0);
% C0 = res*sys.C;
% temp = K*D0;
% B0 = sys.B*(eye(size(temp,1), size(temp,2)) + temp);
% A0 = sys.A + sys.B*(K*C0);