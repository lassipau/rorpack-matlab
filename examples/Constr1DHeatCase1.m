function [x0,Sys,spgrid] = Constr1DHeatCase1(cval,x0fun,N)
% [x0,Sys,spgrid] = Constr1DHeatCase1(x0fun,N)
% 
% Finite Differences approximation of a 1D Heat equation with different
% types of boundary control 
% with Neumann boundary input (at x=0) and
% disturbance and measured temperature output (at x=1)
% Usage: Can also use solely for defining the initial state x0
% cval = thermal diffusivity of the material (assumed constant)
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

spgrid = linspace(0,1,N);
h = 1/(N-1);

ee = ones(N,1);

A = cval*1/h^2*spdiags([ee -2*ee ee],-1:1,N,N);
A(1,2) = cval*2/h^2;
A(N,N-1) = cval*2/h^2;


% Neumann boundary input at x=0, sign is based on the _outwards_ normal derivative
B = sparse([-2/h;zeros(N-1,1)]); 


Bd = sparse([zeros(N-1,1);2/h]); % Neumann boundary disturbance at x=1
C = sparse([zeros(1,N-1), 1]); % Measured temperature at x=1
Cm = sparse([1, zeros(1,N-1)]); % Additional output, measured temperature at x=0

% Case 1 has Neumann boundary conditions at both x=0 and x=1
BCtype = 'NN';

x0 = x0fun(spgrid).';


Sys.A = A;
Sys.B = B;
Sys.Bd = Bd;
Sys.C = C;
Sys.Cm = Cm;
Sys.D = 0;
Sys.Dd = 0;
Sys.Dm = 0;



