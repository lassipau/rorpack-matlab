function [x0,Sys,spgrid,BCtype] = ConstrHeat1DCase4(cfun,x0fun,N)
% [x0,Sys,spgrid] = Constr1DHeatCase2(x0fun,N)
% 
% Finite Differences approximation of a 1D Heat equation with different
% types of distributed or boundary control and observation.
%
% Case 4: Dirichlet boundary control at x=1, regulated output y(t) and a 
% Neumann boundary disturbance at x=0. The system is exponentially stable,
% but does not define a "wellposed" or "regular linear system" on the
% natural state space X=L^2(0,1). Since the current theory does not
% guarantee that the controller designs would work, these are simulations
% are only for experimentation purposes. That is, PROCEED WITH CAUTION! ;)
%
% Usage: Can also use solely for defining the initial state x0
% cfun = spatially variying thermal diffusivity of the material
% x0fun = initial heat profile, function handle
% Sys = system parameters, (sys.A,sys.B,sys.Bd,sys.C,sys.D)
% N = dimension of the approximated system (does not include point x=1)
%
% Cases:
% Case 1: Neumann boundary control at x=0, regulated output y(t) and a 
% Neumann boundary disturbance at x=1. Additional measured output at x=0
%
% Case 2: Neumann boundary control and disturbance at x=0, regulated output 
% y(t) at x=0, Dirichlet boundary condition at x=1
%
% Case 3: Neumann boundary disturbance at x=0, 2 distributed controls and 
% two distributed measurements regulated output y(t) 
%
% Case 4: Dirichlet boundary control at x=1, regulated output y(t) and a 
% Neumann boundary disturbance at x=0. The system is exponentially stable,
% but does not define a "wellposed" or "regular linear system" on the
% natural state space X=L^2(0,1). Since the current theory does not
% guarantee that the controller designs would work, these are simulations
% are only for experimentation purposes. That is, PROCEED WITH CAUTION! ;)
%
% Case 5: Similar to Case 1, but with two boundary inputs and outputs: 
% Neumann boundary control u_1(t) at x=0, and u_2(t) at x=1. Pointwise
% temperature measurements y_1(t) at x=0, and y_2(t) at x=1. Two input
% disturbances w_{dist,1}(t) at x=0 and w_{dist,2}(t) at x=1. The system is
% unstable (eigenvalue at 0), but is impedance passive and can be
% stabilized with negative output feedback.


% The point x=1 has a Dirichlet boundary condition, and is not chosen as a
% variable
spgrid = linspace(0,1,N+1);
h = 1/N;
% Case 4 has a Neumann boundary condition at x=0 and a Dirichlet condition at x=1
BCtype = 'ND';

[A,spgrid] = DiffOp1d(cfun,spgrid,BCtype);

% Dirichlet boundary input at x=1
B = sparse([zeros(N-1,1);cfun(1)/h^2]); 

% A Neumann boundary disturbance at x=0
% This is also affected by the diffusion coefficient at x=0 (according to 
% the approximation in DiffOp1d).
Bd = sparse([2*cfun(0)/h;zeros(N-1,1)]);

C = sparse([1,zeros(1,N-1)]); % Measured temperature at x=0

x0 = x0fun(spgrid(1:N)).';


Sys.A = A;
Sys.B = B;
Sys.Bd = Bd;
Sys.C = C;
Sys.D = 0;
Sys.Dd = 0;