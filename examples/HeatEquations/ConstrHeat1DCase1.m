function [x0,Sys,spgrid,BCtype] = ConstrHeat1DCase1(cfun,x0fun,N)
% Finite Differences approximation of a 1D Heat equation with different
% types of distributed or boundary control and observation.
%
% Case 1: Neumann boundary control at x=0, regulated output y(t) and a 
% Neumann boundary disturbance at x=1. Additional measured output at x=0
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
%   Sys : [struct with fields (at least) A,B,C,D] approximation of the
%   system to be controlled
%
%   spgrid : [1xN matrix] spatial grid of the case
%
%   BCtype : [string] Boundary control type for the case, one of 'NN',
%   'ND', 'DN' or 'DD'.

spgrid = linspace(0,1,N);
h = 1/(N-1);
% Case 1 has Neumann boundary conditions at both x=0 and x=1
BCtype = 'NN';

% modify A and the spatial grid in accordance with the thermal diffusivity
[A,spgrid] = DiffOp1d(cfun,spgrid,BCtype);

% Neumann boundary input at x=0, sign is based on the
% _outwards_ normal derivative
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