function [Sys, x0, phin] = ConstrWave1DCase2(N, Bfun, Bdfun, w0fun, wd0fun)

phin = @(x,n) diag(n>0)*sqrt(2)*sin(n*pi*x);

A = zeros(2*N, 2*N);
B = zeros(2*N, 1);
Bd = zeros(2*N, 1);
C = zeros(1, 2*N);
% Kinf = zeros(1, 2*N);
% Linf = zeros(2*N, 1);
x0 = zeros(2*N, 1);


for n = 1:N
  indran = 2*(n-1)+(1:2);
  
  A(indran,indran) = [0 1;-n^2*pi^2 0];
  B(indran) = [0;integral(@(x) Bfun(x) .* conj(phin(x,n)), 0, 1)];
  
  % Disturbance input
  Bd(indran) = [0;integral(@(x) Bdfun(x) .* conj(phin(x,n)), 0, 1)];
  
  C = B.';
  
  % Compute the approximation of the initial state
  x0(indran) = [integral(@(x) w0fun(x).*conj(phin(x,n)),0,1);...
      integral(@(x) wd0fun(x).*conj(phin(x,n)),0,1)];
end

Sys.A = A;
Sys.B = B;
Sys.Bd = Bd;
Sys.C = C;
Sys.D = 0;
Sys.Dd = 0;
