function [x0,spgrid,sys] = ConstrHeat2DtestEx(cval,x0fun,N)
% [x0,spgrid,sys] = ConstrHeat2Dtest1(c,x0fun,N)
% Usage: Can also use solely for defining the initial state x0
% x0 = initial state
% spgrid = spatial grid of the domain
% sys = system parameters, (sys.A,sys.B,sys.Bd,sys.C,sys.D)

xx = linspace(0,1,N);
yy = xx;
[xxg,yyg] = meshgrid(xx,yy);
spgrid.xx = xxg;
spgrid.yy = yyg;

x0 = reshape(x0fun(xxg,yyg),N^2,1);


if nargout > 2
  h = 1/(N-1);

  ee = ones(N,1);


  Lap1D = 1/h^2*spdiags([ee -2*ee ee],[-1:1],N,N);
  % Neumann boundary conditions: v_0=v_1, v_{N+1}=v_N
  Lap1D(1,1) = -1/h^2; 
  Lap1D(end,end) = -1/h^2;

  A = cval*( kron(Lap1D, speye(N)) + kron(speye(N), Lap1D) );

  % Inputs and outputs

  % Neumann boundary control input, the part of the boundary where x=0
  % Corresponds to delta function on the boundary. Points with y=0 have
  % indices 1,N+1,2N+1,...
  
  
  % Neumann input from the boundary where y=0
  b1 = zeros(N^2,1);
  b1(1+N*(0:(N-1)),1) = 1/h;
  
  B = b1;
  
  % Observation is the integral over the part of the boundary where x=1
  % corresponding indices are N*(N-1)+1..N
  c1 = zeros(1,N^2);
  c1(1,N*(N-1)+(1:N))=h;
  c1(1,N*(N-1)+[1 N])=h/2;
  
  C = c1;
  
  D = 0;
  
  
  % Disturbance signal input operator, Neumann input from the boundary where
  % x=0 (corresponds to "first column" of x in grid form
  %B_d = reshape(2/h*[ones(N,1) zeros(N,N-1)],N^2,1);
  
  % Disturbance signal input operator, Neumann input from the boundary where
  % x=0 and 0<=y<1/2 (corresponds to "first column" of x in grid form
  Bd = reshape(2/h*[[ones(floor(N/2),1); zeros(N-floor(N/2),1)] zeros(N,N-1)],N^2,1);

  % For simplicity: Stabilization with state feedback
  % 1 for interior points, 1/4 for corner points and 1/2 for other boundary
  % points
  K = -cval*pi^2*h^2*ones(1,N^2);
  K([1:N N*(N-1)+(1:N)]) = K([1:N N*(N-1)+(1:N)])/2;
  K([1+N*(0:(N-1)) N+N*(0:(N-1))]) = K([1+N*(0:(N-1)) N+N*(0:(N-1))])/2;
  
  sys.A = A+B*K;
  sys.B = B;
  sys.Bd = Bd;
  sys.C = C;
  sys.D = D;
  
end

