function CLsim = SimCLSys(CLSys,xe0,yref,wdist,tgrid,opts)
% sol = SimCLSys(CLsys,xe0,yref,wdist,tspan,opts)
% Simulate the closed-loop system, version compatible with signals wdist
% and yref, instead of the exosystem as the CL system input
% CLsim = CL state, controlled output, and regulation error, 
% CLsim.xesol = xe(t) evaluated at tgrid
% CLsim.output = y(t) at tgrid
% CLsim.error = e(t) at tgrid
% CLsim.solstruct = solution structure from ode15s for interpolation


% %tspan = [0 26];
% tspan = [0 80];
% opts = odeset('Reltol',1e-6,'Abstol',1e-9);


odefun = @(t,xe) CLSys.Ae*xe+CLSys.Be*[wdist(t);yref(t)];

opts = odeset(opts,'Jacobian',CLSys.Ae);

CLsim.solstruct = ode15s(odefun,[tgrid(1) tgrid(end)],xe0,opts);

xevals = deval(CLsim.solstruct,tgrid);
CLsim.xesol = xevals;
CLsim.output = CLSys.Ce*xevals;
CLsim.error = CLSys.Ce*xevals+CLSys.De*[wdist(tgrid);yref(tgrid)];

