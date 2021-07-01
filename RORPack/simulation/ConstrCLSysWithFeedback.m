function CLSys = ConstrCLSysWithFeedback(ContrSys,Sys,Kappa)

temp = Sys.D*Kappa;
res = inv(eye(size(temp,1), size(temp,2)) - temp);
D0 = res*Sys.D;
Dd0 = res*Sys.Dd;
Bd0 = Sys.Bd + Sys.B*(Kappa*Dd0);
C0 = res*Sys.C;
temp = Kappa*D0;
B0 = Sys.B*(eye(size(temp,1), size(temp,2)) + temp);
A0 = Sys.A + Sys.B*(Kappa*C0);

CLSys.Ae = [A0, B0*ContrSys.K;ContrSys.G2*C0, ContrSys.G1 + ContrSys.G2*(D0*ContrSys.K)];
CLSys.Be = [Bd0, -(B0*Kappa);ContrSys.G2*Dd0, -ContrSys.G2*res];
CLSys.Ce = [C0, D0*ContrSys.K];
CLSys.De = [Dd0, -res];

% CLSys.Ae = [Sys.A Sys.B*ContrSys.K;ContrSys.G2*Sys.C ContrSys.G1+ContrSys.G2*Sys.D*ContrSys.K];
% CLSys.Be = [Sys.Bd zeros(size(Sys.Bd,1),size(Sys.C,1));zeros(size(ContrSys.G1,1),size(Sys.Bd,2)) -ContrSys.G2];
% CLSys.Ce = [Sys.C Sys.D*ContrSys.K];
% CLSys.De = [zeros(size(Sys.C,1),size(Sys.Bd,2)) -eye(size(Sys.C,1))];