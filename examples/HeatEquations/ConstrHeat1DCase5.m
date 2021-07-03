function [x0,Sys,spgrid,BCtype] = ConstrHeat1DCase5(cfun,x0fun,N)
% [x0,Sys,spgrid] = Constr1DHeatCase5(x0fun,N)
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
% Case 5 has Neumann boundary conditions at both x=0 and x=1
BCtype = 'NN';

[A,spgrid] = DiffOp1d(cfun,spgrid,BCtype);


% Neumann boundary input at x=0 (input u_1(t)) and x=1 (input u_2(t)) 
% Signs are based on the _outwards_ normal derivatives
B = [[2*cfun(0)/h;zeros(N-1,1)],[zeros(N-1,1);-2*cfun(1)/h]]; 

% Two input disturbances
Bd = B;

% Measured temperatures at x=0 (y_1(t)) and x=1 (y_2(t))
C = sparse([1, zeros(1,N-1);zeros(1,N-1), 1]); 


x0 = x0fun(spgrid).';

Sys.A = A;
Sys.B = B;
Sys.Bd = Bd;
Sys.C = C;
Sys.D = zeros(2,2);
Sys.Dd = zeros(2,1);




