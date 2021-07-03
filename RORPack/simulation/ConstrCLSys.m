function CLSys = ConstrCLSys(Sys,ContrSys)
% function CLsys = ConstrCLSys(Sys,ContrSys)
%
% Construct the closed-loop system of the system 'Sys' and the dynamic
% error feedback controller 'ContrSys'. The inputs of the closed-loop is
% of the form [wdist(t);yref(t)] where wdist(t) is the disturbance signal 
% and yref(t) is the reference signal. The use of the consistency check 
% routine "SysConsistent" is recommended before calling ConstrCLSys.
%
% Inputs: 
% 'Sys' = [struct] The linear system with parameters (Sys.A,Sys.B,Sys.Bd,Sys.C,Sys.D,Sys.Dd)
% 'ContrSys' = [struct] The controller with parameters (ContrSys.G1,ContrSys.G2,ContrSys.K,ContrSys.Dc)
%
% Outputs:
% 'CLSys' = [struct] The closed-loop system with parameters (CLsys.Ae,CLsys.Be,CLsys.Ce,CLsys.De) 

if ~isfield(Sys,'Bd')
    error('Disturbance input operator "Bd" is not defined! (Use "SysConsistent" before closed-loop construction).')
elseif ~isfield(Sys,'Dd')
    Sys.Dd = zeros(size(Sys.C,1),size(Sys.Bd,2));
end


dimX = size(Sys.A,1);
dimU = size(Sys.B,2);
dimUd = size(Sys.Bd,2);
dimY = size(Sys.C,1);


% Version with possible controller feedthrough term D_c

if ~isfield(ContrSys,'Dc') || isequal(ContrSys.Dc,zeros(dimU,dimY))
    CLSys.Ae = [Sys.A, Sys.B*ContrSys.K;ContrSys.G2*Sys.C, ContrSys.G1+ContrSys.G2*Sys.D*ContrSys.K];
    CLSys.Be = [Sys.Bd, zeros(dimX,dimY);ContrSys.G2*Sys.Dd, -ContrSys.G2];
    CLSys.Ce = [Sys.C, Sys.D*ContrSys.K];
    CLSys.De = [Sys.Dd, -eye(dimY)];
    
    % An additional output operator for computing the control input u(t) based
    % on the closed-loop state x_e(t)
    CLSys.Ke = [zeros(size(ContrSys.K,1),size(Sys.A,1)), ContrSys.K];
    CLSys.Due = zeros(dimU,dimUd+dimY);
else
    
    % Construct the closed-loop system by interpreting the controller 
    %feedthrough as a preliminary output feedback to the system 'Sys'.
    SysFB = SysOutputFeedback(Sys, ContrSys.Dc);
    
    CLSys.Ae = [SysFB.A, SysFB.B*ContrSys.K;ContrSys.G2*SysFB.C, ContrSys.G1+ContrSys.G2*SysFB.D*ContrSys.K];
    CLSys.Be = [SysFB.Bd, -SysFB.B*ContrSys.Dc; ContrSys.G2*SysFB.Dd, -ContrSys.G2/(eye(dimY)-Sys.D*ContrSys.Dc)];
        
    CLSys.Ce = [SysFB.C, SysFB.D*ContrSys.K];
    CLSys.De = [SysFB.Dd, -inv(eye(dimY)-Sys.D*ContrSys.Dc)];
    
    % An additional output operator for computing the control input u(t) based
    % on the closed-loop state x_e(t)
    CLSys.Ke = [ContrSys.Dc*SysFB.C, (eye(dimU)-ContrSys.Dc*Sys.D)\ContrSys.K];
    CLSys.Due = [ContrSys.Dc*SysFB.Dd,-ContrSys.Dc/(eye(dimY)-Sys.D*ContrSys.Dc)];
end