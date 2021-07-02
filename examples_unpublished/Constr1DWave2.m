function [x0,Sys,phin] = Constr1DWave(w0fun,wd0fun,N)
% [x0,Sys,phin] = Constr1DWave(x0fun,N)
% 
% Modal approximation of a 1D Wave equation with Neumann boundary input (at x=1) and
% disturbance (at x=0), and Dirichlet outputs y(t)=(w(1,t),w_t(1,t))^T
% Usage: Can also use solely for defining the initial state x0
% w0fun = initial displacement, function handle
% wd0fun = initial velocity, function handle
% sys = system parameters, (sys.A,sys.B,sys.Bd,sys.C,sys.D)
% phin = Normalized eigenfunction of the Laplacian

% N = 14;

phin = @(x,n) (n==0)*ones(size(x)) + diag(n>0)*sqrt(2)*cos(n*pi*x);

A = zeros(2*N);
B = zeros(2*N,1);
Bd = zeros(2*N,1);
C = zeros(1,2*N);
C2 = zeros(1,2*N);

x0 = zeros(2*N,1);

for n_val = 0:(N-1)
  indran = 2*n_val+(1:2);
  A(indran,indran) = [0 1;-n_val^2*pi^2 0];
  B(indran) = [0;phin(1,n_val)];
  Bd(indran) = [0;phin(0,n_val)];
%   C(indran) = [phin(1,n_val) 0];
%   C2(indran) = [phin(0,n_val) 0]; % Displacement at x=0

  C(indran) = [0 phin(0,n_val)]; % Velocity at x=0
%   C2(indran) = [0 phin(1,n_val)];

  x0(indran) = [integral(@(x) w0fun(x).*conj(phin(x,n_val)),0,1);integral(@(x) wd0fun(x).*conj(phin(x,n_val)),0,1)];
end

C2 = [1, 0, zeros(1,2*(N-1))]; % Average displacement

Sys.A = sparse(A);

Sys.B = B;
Sys.Bd = Bd;
% Sys.C = [C;C2];
Sys.C = C;
Sys.Cm = C2;
% Sys.D = [0;0];
Sys.D = 0;
% Sys.Dm = 0;

% Sys.A = A-1/2*B*C2;

%  Sys.A = A-1/2*B*C2;
% % Sys.B = B;
% %Sys.A = A-1/2*B*(C+C2);
% Sys.C = C;
% Sys.C2 = C2;
% Sys.D = 0;



