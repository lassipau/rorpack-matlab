function [x0,Sys,phin,Kinf,Linf] = ConstrWave1DCase1(w0fun,wd0fun,N) 
% Modal approximation of a 1D Wave equation with Neumann boundary input
% (at x=1) and disturbance (at x=0), and Dirichlet velocity 
% output y(t)=w_t(0,t)^T. The additional measured output 
% y_m(t) = \int_0^1 y(xi,t) dxi can be used in stabilising the zero 
% eigenvalue of the system.
%
% Usage: Can also be used solely for defining the initial state x0
% w0fun = initial displacement, function handle
% wd0fun = initial velocity, function handle
% sys = system parameters, (sys.A,sys.B,sys.Bd,sys.C,sys.D)
% phin = Normalized eigenfunction of the Laplacian
%
% Also designs the "infinite parts" of the exponentially stabilizing
% feedback gains L and K

phin = @(x,n) (n==0)*ones(size(x)) + diag(n>0)*sqrt(2)*cos(n*pi*x);

A = zeros(2*N);
B = zeros(2*N,1);
Bd = zeros(2*N,1);
C = zeros(1,2*N);

Kinf = zeros(1,2*N);
Linf = zeros(2*N,1);

% C2 = zeros(1,2*N);

x0 = zeros(2*N,1);

for n_val = 0:(N-1)
  indran = 2*n_val+(1:2);
  A(indran,indran) = [0 1;-n_val^2*pi^2 0];
  B(indran) = [0;phin(1,n_val)];
  Bd(indran) = [0;phin(0,n_val)];
%   C(indran) = [phin(1,n_val) 0];
%   C(indran) = [phin(0,n_val) 0];
  C(indran) = [0, phin(0,n_val)];
%   C2(indran) = [0 phin(1,n_val)];

  Kinf(indran) = [0, phin(1,n_val)];
  Linf(indran) = [0;phin(0,n_val)];

  % Compute the approximation of the initial state
  x0(indran) = [integral(@(x) w0fun(x).*conj(phin(x,n_val)),0,1);...
      integral(@(x) wd0fun(x).*conj(phin(x,n_val)),0,1)];
end

Sys.A = sparse(A);
Sys.B = B;
Sys.Bd = Bd;
Sys.C = C;
% Sys.Cm = C2;
Sys.D = 0;
% Sys.Dm = 0;
Sys.Dd = zeros(size(Sys.C,1),size(Sys.Bd,2));


% Cm = measure the average displacement on the wave profile, used in
% prestabilization of the eigenvalue \lambda=0
Sys.Cm = [[1,0],zeros(1,2*(N-1))];
Sys.Dm = 0;

