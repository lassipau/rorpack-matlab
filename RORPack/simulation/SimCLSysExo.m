function CLsim = SimCLSysExo(CLsys,xe0,exosys,v0,tgrid,opts)
% currently unused in the given example files
% sol = SimCLSys(CLsys,xe0,exosys,v0,tspan,opts)
% Simulate the closed-loop system
% sol = CL+exo state and regulation error, CL.xasol, CLsim.error

Aaug = [CLsys.Ae CLsys.Be;zeros(length(v0),length(x0)) exosys.S];
x0aug = [xe0;v0];


% %tspan = [0 26];
% tspan = [0 80];
% opts = odeset('Reltol',1e-6,'Abstol',1e-9);

odefun = @(t,xa) Aaug*xa;

sol = ode15s(odefun,[tgrid(1) tgrid(2)],x0aug,opts);

xavals = deval(sol,tgrid);
CLsim.xesol = xavals(1:length(xe0),:);
CLsim.vsol = xavals((length(xe0)+1):end,:);
CLsim.error = [CLsys.Ce CLsys.De]*xavals;
