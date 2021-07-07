function [x0,Sys,spgrid,BCtype] = ConstrHeat1DCase3(cfun,x0fun,N,IB1,IB2,IC1,IC2)
% Finite Differences approximation of a 1D Heat equation with different
% types of distributed or boundary control and observation.
%
% Case 3: Neumann boundary disturbance at x=0, two distributed controls and 
% two distributed outputs y(t). The controls affect the intervals 
% 'IB1' = [a_1,b_1] and 'IB2' = [a_2,b_2], and the measurements are
% the averages of the temperatures on the intervals 'IC1' = [c_1,d_1] and 
% 'IC2' = [c_2,d_2]. If these parameters are not given, then default 
% configuration IB1 = [.3, .4], IB2 = [.6, .7], IC1 = [.1, .2], and 
% IC2 = [.8, .9] is used.
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
%   IB1,IB2,IC1,IC2 : [1x2 matrix] control and measurement intervals
%   (optional)
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
% disturbances w_{dist,1}(t) at x=0 and w_{dist,2}(t) at x=1, and a third
% distributed disturbance with profile "Bd_profile" (function). The system 
% is unstable (eigenvalue at 0), but is impedance passive and can be
% stabilized with negative output feedback.



if nargin <= 3
  IB1 = [.3, .4];
  IB2 = [.6, .7];
  IC1 = [.1, .2];
  IC2 = [.8, .9];
end

% Check that the intervals are defined correctly
if ~isempty(find([IB1(:);IB2(:);IC1(:);IC2(:)]<0))...
        || ~isempty(find([IB1(:);IB2(:);IC1(:);IC2(:)]>1))
  error(...
      'All limits of the intervals IB1, IB2, IC1, and IC2 need to be on [0,1].')
elseif IB1(2)<=IB1(1) || IB2(2)<=IB2(1) ||...
        IC1(2)<=IC1(1) || IC2(2)<=IC2(1)
  error(...
      'All intervals IB1, IB2, IC1, and IC2 need to be of positive lengths.')
end

% The point x=1 has a Dirichlet boundary condition
% and is not chosen as a variable
spgrid = linspace(0,1,N+1);
h = 1/N;
% Case 3 has a Neumann boundary condition at x=0
% and a Dirichlet condition at x=1
BCtype = 'ND';

% modify A and the spatial grid in accordance with the thermal diffusivity
[A,spgrid] = DiffOp1d(cfun,spgrid,BCtype);

% Two distributed inputs on the intervals IB1 and IB2
B1 = 1/(IB1(2)-IB1(1))*(spgrid(1:N)>=IB1(1) & spgrid(1:N)<=IB1(2));
B1 = B1(:);
B2 = 1/(IB2(2)-IB2(1))*(spgrid(1:N)>=IB2(1) & spgrid(1:N)<=IB2(2));
B2 = B2(:);

% Neumann boundary disturbance at x=0, sign is based on the
% _outwards_ normal derivative
Bd = sparse([2*cfun(0)/h;zeros(N-1,1)]); 

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