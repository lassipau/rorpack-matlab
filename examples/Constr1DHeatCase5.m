function [x0,Sys,spgrid,BCtype] = Constr1DHeatCase5(x0fun,N)
% [x0,Sys,spgrid] = Constr1DHeatCase5(x0fun,N)
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

% Petteri's parameters: cval = 1, gamma = pi^2+1
cval = 1;
gamma = pi^2+1;

A = cval*1/h^2*spdiags([ee -2*ee ee],-1:1,N,N);
A(1,2) = cval*2/h^2;
A(N,N-1) = cval*2/h^2;
A = A + gamma*spdiags(ee,0,N,N);


% Neumann boundary input at x=0 (input u_1(t)) and x=1 (input u_2(t)) 
% Signs are based on the _outwards_ normal derivatives
B = [[-2/h;zeros(N-1,1)],[zeros(N-1,1);2/h]]; 


IC1 = [.0, .25];
IC2 = [.5, .75];

% Two distributed outputs on the intervals IC1 and IC2
weights = [1/2,ones(1,N-2),1/2];
C1 = h/(IC1(2)-IC1(1))*(spgrid>=IC1(1) & spgrid<=IC1(2)).*weights;
C2 = h/(IC2(2)-IC2(1))*(spgrid>=IC2(1) & spgrid<=IC2(2)).*weights;

C = [C1;C2]; 

% Case 5 has Neumann boundary conditions at both x=0 and x=1
BCtype = 'NN';

x0 = x0fun(spgrid).';


Sys.A = A;
Sys.B = B;
Sys.Bd = zeros(N,1);
Sys.C = C;
Sys.D = 0;
Sys.Dd = zeros(2,1);




