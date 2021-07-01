function stabmarg = CLStabMargin(CLSys,raise_error)
% stabmarg = CLStabMargin(CLSys)
%
% Return the stability margin of the CL system
% Negative stability margin means the CL system is unstable
% 
% If 'raise_error' is 'true', then an unstable closed-loop system causes an
% error. By default, 'raise_error' is 'false', and an unstable closed-loop
% system causes a warning.

stabmarg = -max(real(eig(full(CLSys.Ae))));
if stabmarg < 0
    if nargin==2 && raise_error
        error("Stability margin < 0. The closed-loop system is unstable!")
    else
        warning("Stability margin < 0. The closed-loop system is unstable!")
    end
end