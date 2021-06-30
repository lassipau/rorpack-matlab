function [x0,spgrid,sys] = ConstrHeat2Dtest2(N,M,x0fun,cval)
% [x0,spgrid,sys] = ConstrHeat2Dtest1(c,x0fun,N)
% Usage: Can also use solely for defining the initial state x0
% x0 = initial state
% spgrid = spatial grid of the domain
% sys = system parameters, (sys.A,sys.B,sys.Bd,sys.C,sys.D)

xx = linspace(0,1,N);
yy = linspace(0,1,M);
[xxg,yyg] = meshgrid(xx,yy);
spgrid.xx = xxg;
spgrid.yy = yyg;

x0 = reshape(x0fun(xxg,yyg),N*M,1);

if nargout > 2
  dx = 1/(N-1);
  dy = 1/(M-1);

  ex = ones(N,1);
  ey = ones(M,1);

  Dxx = spdiags([ex -2*ex ex],[-1:1],N,N);
  Dyy = spdiags([ey -2*ey ey],[-1:1],M,M);
  
  % Neumann boundary conditions: v_0=v_1, v_{N+1}=v_N
  Dxx(1,1) = -1; 
  Dxx(end,end) = -1;
  Dyy(1,1) = -1;
  Dyy(end,end) = -1;
  
  A = cval * full(kron(Dyy, speye(N)) + kron(speye(M), Dxx))/(dx*dy);
  % Inputs and outputs

  % Neumann boundary control input, the part of the boundary where x=0
  % Corresponds to delta function on the boundary. Points with y=0 have
  % indices 1,N+1,2N+1,...
  
  
  % Neumann input from the boundary where x=0
  b1 = zeros(N*M,1);
  b1(1:N:N*M,1) = 1/dy;
  
  B = b1;
  
  % Observation is the integral over the part of the boundary where y=1
  % corresponding indices are N*(N-1)+1..N
  c1 = zeros(1,N*M);
  c1(1,(N-1)*M+1:end) = dx;
  c1(1,(N-1)*M+1) = dx / 2; 
  c1(1, end) = dx / 2;
  
  C = c1;
  
  D = 0;
  
  
  % Disturbance signal input operator, Neumann input from the boundary where
  % x=0 (corresponds to "first column" of x in grid form
  %B_d = reshape(2/h*[ones(N,1) zeros(N,N-1)],N^2,1);
  
  % Disturbance signal input operator, Neumann input from the boundary where
  % x=0 and 0<=y<1/2 (corresponds to "first column" of x in grid form
  Bd = reshape(2/dx*[[ones(floor(N/2),1); zeros(N-floor(N/2),1)] zeros(N,M-1)],N*M,1);

  sys.A = A;
  sys.B = B;
  sys.Bd = Bd;
  sys.C = C;
  sys.D = D;
  
end

