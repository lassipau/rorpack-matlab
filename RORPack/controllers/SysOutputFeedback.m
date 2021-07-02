function SysFB = SysOutputFeedback(Sys, K)
% Applies output feedback y(t)=K*u(t)+u_new(t) to the system 'Sys'
% (A+BK(I+DK)^(-1)C,B(I+KD)^(-1),(I+DK)^(-1)C,(I+DK)^(-1)D) 
% obtained from (A,B,C,D) with output feedback u(t)=Ky(t)
%
% Input parameters:
% Sys = The linear system (A,B,Bd,C,D,Dd).
% K = The output feedback matrix.
%
% Output parameters
% SysFB = The system after output feedback.

dimU = size(Sys.D,2);
dimY = size(Sys.D,1);

if ~isequal(size(K),[dimU,dimY])
    error('Output feedback cannot be applied, the matrix K has incorrect dimensions!')
elseif abs(det(eye(dimY)-Sys.D*K))<1e-10
    error('Output feedback cannot be applied, the matrix I-DK is singular! (with tolerance 1e-10)')
end

% Compute the modified main system operators 
SysFB.B = Sys.B/(eye(dimU) - K*Sys.D);
SysFB.A = Sys.A + SysFB.B*(K*Sys.C);
SysFB.C = (eye(dimY) - Sys.D*K)\Sys.C;
SysFB.D = (eye(dimY) - Sys.D*K)\Sys.D;

% Compute the modified disturbance input and feedthrough operators.
% If only one of Bd or Dd is defined, also define the other one zero (will 
% be used later).
if isfield(Sys,'Dd')
    SysFB.Dd = (eye(dimY) - Sys.D*K)\Sys.Dd;
    if isfield(Sys,'Bd')
        SysFB.Bd = Sys.Bd + SysFB.B*K*Sys.Dd;
    else 
        SysFB.Bd = zeros(size(Sys.A,1),size(Sys.Dd,2));
        Sys.Bd = zeros(size(Sys.A,1),size(Sys.Dd,2));
    end
elseif isfield(Sys,'Bd')
    SysFB.Bd = Sys.Bd;
    SysFB.Dd = zeros(dimY,size(Sys.Bd,2));
    Sys.Dd = zeros(dimY,size(Sys.Bd,2));
end


% Compute the modified measured output and feedthrough operators
% (only done if 'Cm' is defined).
if isfield(Sys,'Cm')
    if ~isfield(Sys,'Dm')
        Sys.Dm = zeros(size(Cm,1),dimU);
    end
    SysFB.Cm = Sys.Cm+Sys.Dm*K*SysFB.C;
    SysFB.Dm = Sys.Dm+Sys.Dm*K*SysFB.D;
    
    % If the system has a disturbance input, also compute modified 'Dmd'.
    if isfield(Sys,'Dd')
        if ~isfield(Sys,'Dmd')
            Sys.Dmd = zeros(size(Sys.Cm,1),size(Sys.Dd,2));
        end
        SysFB.Dmd = Sys.Dmd+Sys.Dm*K*SysFB.Dd;
    end
end
