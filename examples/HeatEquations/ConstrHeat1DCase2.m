function [x0,Sys,spgrid,BCtype] = ConstrHeat1DCase2(cfun,x0fun,N)
% Finite Differences approximation of a 1D Heat equation with different
% types of boundary control 
% with Neumann boundary input (at x=0) and
% disturbance and measured temperature output (at x=1)
%
% Case 2: Neumann boundary control and disturbance at x=0, regulated output 
% y(t) at x=0, Dirichlet boundary condition at x=1
%
% Usage: Can also use solely for defining the initial state x0
% Inputs:
%   cfun : [1x1 function_handle] spatially varying
%   thermal diffusivity of the material
%
%   x0fun : [1x1 function_handle] initial heat profile
%
%   N : [integer] Dimension of the approximated system
%
% Outputs:
%   x0 : [Nx1 matrix] initial state of the plant
%
%   Sys : [struct with fields A,B,C,D] approximation of the
%   system to be controlled
%
%   spgrid : [1xN matrix] spatial grid of the case
%
%   BCtype : [string] Boundary control type for the case, one of 'NN',
%   'ND', 'DN' or 'DD'.

% The point x=1 has a Dirichlet boundary condition, and is not chosen as a
% variable
spgrid = linspace(0,1,N+1);
h = 1/N;
% Case 2 has a Neumann boundary condition at x=0
% and a Dirichlet condition at x=1
BCtype = 'ND';

% modify A and the spatial grid in accordance with the thermal diffusivity
[A,spgrid] = DiffOp1d(cfun,spgrid,BCtype);

% Neumann boundary input at x=0, sign is based on the
% _outwards_ normal derivative
% This is also affected by the diffusion coefficient at x=0 (according to 
% the approximation in DiffOp1d).
B = sparse([2*cfun(0)/h;zeros(N-1,1)]);

Bd = B; % Neumann boundary disturbance at x=0

C = sparse([1,zeros(1,N-1)]); % Measured temperature at x=0

x0 = x0fun(spgrid(1:N)).';


Sys.A = A;
Sys.B = B;
Sys.C = C;
Sys.D = 0;
Sys.Bd = Bd;
Sys.Dd = 0;