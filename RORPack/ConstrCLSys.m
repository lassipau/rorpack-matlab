function CLSys = ConstrCLSys(Sys,ContrSys)
% CLsys = ConstrCLSys(Sys,ContrSys)
%
% Construct the closed-loop system
% Version corresponding to the inputs wdist(t) and yref(t) instead of the
% exosystem state
% CLsys = CL system parameters, (CLsys.Ae,CLsys.Be,CLsys.Ce,CLsys.De) 

CLSys.Ae = [Sys.A Sys.B*ContrSys.K;ContrSys.G2*Sys.C ContrSys.G1+ContrSys.G2*Sys.D*ContrSys.K];
CLSys.Be = [Sys.Bd zeros(size(Sys.Bd,1),size(Sys.C,1));zeros(size(ContrSys.G1,1),size(Sys.Bd,2)) -ContrSys.G2];
CLSys.Ce = [Sys.C Sys.D*ContrSys.K];
CLSys.De = [zeros(size(Sys.C,1),size(Sys.Bd,2)) -eye(size(Sys.C,1))];