function [x0,Sys,spgrid,BCtype] = ConstrHeat1DCase1(cfun,x0fun,N)
% [x0,Sys,spgrid] = Constr1DHeatCase1(x0fun,N)
% 
% Finite Differences approximation of a 1D Heat equation with different
% types of boundary control 
% with Neumann boundary input (at x=0) and
% disturbance and measured temperature output (at x=1)
% Usage: Can also use solely for defining the initial state x0
% cfun = spatially variying thermal diffusivity of the material
% x0fun = initial heat profile, function handle
% Sys = system parameters, (sys.A,sys.B,sys.Bd,sys.C,sys.D)
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


spgrid = linspace(0,1,N);
h = 1/(N-1);
% Case 1 has Neumann boundary conditions at both x=0 and x=1
BCtype = 'NN';

[A,spgrid] = DiffOp1d(cfun,spgrid,BCtype);

% Neumann boundary input at x=0, sign is based on the _outwards_ normal derivative
% This is also affected by the diffusion coefficient at x=0 (according to 
% the approximation in DiffOp1d).
B = sparse([2*cfun(0)/h;zeros(N-1,1)]);


Bd = sparse([zeros(N-1,1);2*cfun(1)/h]); % Neumann boundary disturbance at x=1
C = sparse([zeros(1,N-1), 1]); % Measured temperature at x=1
Cm = sparse([1, zeros(1,N-1)]); % Additional output, measured temperature at x=0

x0 = x0fun(spgrid).';


Sys.A = A;
Sys.B = B;
Sys.Bd = Bd;
Sys.C = C;
Sys.Cm = Cm;
Sys.D = 0;
Sys.Dd = 0;
Sys.Dm = 0;