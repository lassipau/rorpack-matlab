function stabmarg = CLStabMargin(CLSys)
% stabmarg = CLStabMargin(CLSys)
%
% Return the stability margin of the CL system
% Negative stability margin means the CL system is unstable

stabmarg = -max(real(eig(full(CLSys.Ae))));
if stabmarg < 0
    warning("Stability margin < 0 indicates the system is unstable.")
end