function CLsim = SimCLSys(CLSys,xe0,yref,wdist,tgrid,opts)
% sol = SimCLSys(CLsys,xe0,yref,wdist,tspan,opts)
% Simulate the closed-loop system
% Parameters:
% CLSys = the closed-loop system (structure)
% xe0 = the initial state of the closed-loop system
% yref = the reference signal
% wdist = the disturbance signal
% tgrid = the time points where the solution should be evaluated
% opts = additional options to be passed to the ode15s solver (the Jacobian
%        is defined by default in the code). For higher accuracy of the
%        solver, one can set:
%        opts = odeset('Reltol',1e-6,'Abstol',1e-9);
%
% Output CLsim is a structure with the following fields: 
% CLsim = CL state, controlled output, and regulation error, 
% CLsim.xesol = closed-loop state xe(t) evaluated at tgrid
% CLsim.output = the measured output y(t) at tgrid
% CLsim.error = the tracking error e(t)=y(t)-yref(t) at tgrid
% CLsim.control = the control input u(t) at tgrid
% CLsim.solstruct = solution structure from ode15s for interpolation

odefun = @(t,xe) CLSys.Ae*xe+CLSys.Be*[wdist(t);yref(t)];
opts = odeset(opts,'Jacobian',CLSys.Ae);

CLsim.solstruct = ode15s(odefun,[tgrid(1) tgrid(end)],xe0,opts);

xevals = deval(CLsim.solstruct,tgrid);
CLsim.xesol = xevals;
CLsim.output = CLSys.Ce*xevals;
CLsim.error = CLSys.Ce*xevals+CLSys.De*[wdist(tgrid);yref(tgrid)];
CLsim.control = CLSys.Ke*xevals;
