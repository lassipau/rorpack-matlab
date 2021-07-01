function CLSys = ConstrCLSys(Sys,ContrSys)
% CLsys = ConstrCLSys(Sys,ContrSys)
%
% Construct the closed-loop system
% The inputs of the closed-loop system are the disturbance signal wdist(t) 
% and the reference signal yref(t)
% CLsys = CL system parameters, (CLsys.Ae,CLsys.Be,CLsys.Ce,CLsys.De) 

CLSys.Ae = [Sys.A, Sys.B*ContrSys.K;ContrSys.G2*Sys.C, ContrSys.G1+ContrSys.G2*Sys.D*ContrSys.K];
CLSys.Be = [Sys.Bd, zeros(size(Sys.Bd,1),size(Sys.C,1));zeros(size(ContrSys.G1,1),size(Sys.Bd,2)), -ContrSys.G2];
CLSys.Ce = [Sys.C, Sys.D*ContrSys.K];
CLSys.De = [zeros(size(Sys.C,1),size(Sys.Bd,2)), -eye(size(Sys.C,1))];

% An additional output operator for computing the control input u(t) based
% on the closed-loop state x_e(t)
CLSys.Ke = [zeros(size(ContrSys.K,1),size(Sys.A,1)), ContrSys.K];
