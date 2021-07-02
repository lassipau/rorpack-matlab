function CLSys = ConstrCLSys(Sys,ContrSys)
% CLsys = ConstrCLSys(Sys,ContrSys)
%
% Construct the closed-loop system
% The inputs of the closed-loop system are the disturbance signal wdist(t) 
% and the reference signal yref(t)
% CLsys = CL system parameters, (CLsys.Ae,CLsys.Be,CLsys.Ce,CLsys.De) 

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