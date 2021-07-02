function [x0,Sys,spgrid,BCtype] = Constr1DHeatCase5(cfun,x0fun,N)
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

spgrid = linspace(0,1,N);
h = 1/(N-1);
% Case 5 has Neumann boundary conditions at both x=0 and x=1
BCtype = 'NN';

[A,spgrid] = DiffOp1d(cfun,spgrid,BCtype);

% Petteri's parameters: cval = 1, gamma = pi^2+1
gamma = pi^2+1;

A = A + gamma*spdiags(ones(N,1),0,N,N);


% Neumann boundary input at x=0 (input u_1(t)) and x=1 (input u_2(t)) 
% Signs are based on the _outwards_ normal derivatives
B = [[2*cfun(0)/h;zeros(N-1,1)],[zeros(N-1,1);-2*cfun(1)/h]]; 


IC1 = [.0, .25];
IC2 = [.5, .75];

% Two distributed outputs on the intervals IC1 and IC2
weights = [1/2,ones(1,N-2),1/2];
C1 = h/(IC1(2)-IC1(1))*(spgrid>=IC1(1) & spgrid<=IC1(2)).*weights;
C2 = h/(IC2(2)-IC2(1))*(spgrid>=IC2(1) & spgrid<=IC2(2)).*weights;

C = [C1;C2];

x0 = x0fun(spgrid).';


Sys.A = A;
Sys.B = B;
Sys.Bd = zeros(N,1);
Sys.C = C;
Sys.D = zeros(2,2);
Sys.Dd = zeros(2,1);




