function [x0,Sys,spgrid,BCtype] = Constr1DHeatCase3(cfun,x0fun,N,IB1,IB2,IC1,IC2)
% [x0,Sys,spgrid] = Constr1DHeatCase3(x0fun,N)
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
% two distributed measurements regulated output y(t). The controls affect
% the intervals 'IB1' = [a_1,b_1] and 'IB2' = [a_2,b_2], and the
% measurements are the averages of the temperatures on the intervals 
% 'IC1' = [c_1,d_1] and 'IC2' = [c_2,d_2]. If these parameters are not
% given, then default configuration   IB1 = [.3, .4], IB2 = [.6, .7], 
% IC1 = [.1, .2], and IC2 = [.8, .9] is used.

if nargin <= 3
  IB1 = [.3, .4];
  IB2 = [.6, .7];
  IC1 = [.1, .2];
  IC2 = [.8, .9];
end

% Check that the intervals are defined correctly
if ~isempty(find([IB1(:);IB2(:);IC1(:);IC2(:)]<0)) || ~isempty(find([IB1(:);IB2(:);IC1(:);IC2(:)]>1))
  error('All limits of the intervals IB1, IB2, IC1, and IC2 need to be on [0,1].')
elseif IB1(2)<=IB1(1) || IB2(2)<=IB2(1) || IC1(2)<=IC1(1) || IC2(2)<=IC2(1)
  error('All intervals IB1, IB2, IC1, and IC2 need to be of positive lengths.')
end

% The point x=1 has a Dirichlet boundary condition, and is not chosen as a
% variable
spgrid = linspace(0,1,N+1);
h = 1/N;
% Case 3 has a Neumann boundary condition at x=0 and a Dirichlet condition at x=1
BCtype = 'ND';

[A,spgrid] = DiffOp1d(cfun,spgrid,BCtype);

% Two distributed inputs on the intervals IB1 and IB2
B1 = 1/(IB1(2)-IB1(1))*(spgrid(1:N)>=IB1(1) & spgrid(1:N)<=IB1(2));
B1 = B1(:);
B2 = 1/(IB2(2)-IB2(1))*(spgrid(1:N)>=IB2(1) & spgrid(1:N)<=IB2(2));
B2 = B2(:);

% Neumann boundary disturbance at x=0, sign is based on the _outwards_ normal derivative
Bd = sparse([-2/h;zeros(N-1,1)]); 

% Two distributed outputs on the intervals IC1 and IC2
C1 = h/(IC1(2)-IC1(1))*(spgrid(1:N)>=IC1(1) & spgrid(1:N)<=IC1(2));
C2 = h/(IC2(2)-IC2(1))*(spgrid(1:N)>=IC2(1) & spgrid(1:N)<=IC2(2));

x0 = x0fun(spgrid(1:N)).';

Sys.A = A;
Sys.B = [B1,B2];
Sys.Bd = Bd;
Sys.C = [C1;C2];
Sys.D = zeros(2,2);
Sys.Dd = zeros(2,1);