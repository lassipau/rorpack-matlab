function [x0,Sys,spgrid,BCtype] = Constr1DHeatCase2(cfun,x0fun,N)
% [x0,Sys,spgrid] = Constr1DHeatCase2(x0fun,N)
% 
% Finite Differences approximation of a 1D Heat equation with different
% types of boundary control 
% with Neumann boundary input (at x=0) and
% disturbance and measured temperature output (at x=1)
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

% The point x=1 has a Dirichlet boundary condition, and is not chosen as a
% variable
spgrid = linspace(0,1,N+1);
h = 1/N;
% Case 2 has a Neumann boundary condition at x=0 and a Dirichlet condition at x=1
BCtype = 'ND';

[A,spgrid] = DiffOp1d(cfun,spgrid,BCtype);

% Neumann boundary input at x=0, sign is based on the _outwards_ normal derivative
B = sparse([2/h;zeros(N-1,1)]);

Bd = B; % Neumann boundary disturbance at x=0

C = sparse([1,zeros(1,N-1)]); % Measured temperature at x=0

x0 = x0fun(spgrid(1:N)).';


Sys.A = A;
Sys.B = B;
Sys.C = C;
Sys.D = 0;
Sys.Bd = Bd;
Sys.Dd = 0;