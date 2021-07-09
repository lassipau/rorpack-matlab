function [x0,spgrid,Sys] = ConstrTimoshenkoLHMNC18(w0fun,wd0fun,phi0fun,phid0fun,N)
% function [x0,spgrid,Sys] = ConstrTimoshenkoLHMNC18(x10fun,x20fun,x30fun,x40fun,N)
% [x0,spgrid,sys] = ConstrLHMNC18(x10fun,x20fun,x30fun,x40fun,N)
%
% Construct the Timoshenko beam model from the conference article by
% Paunonen, Le Gorrec and Ramirez at LHMNC in 2018.
% Usage: Can also use solely for defining the initial state x0
% x0 = initial state
% spgrid = spatial grid of the domain
% sys = system parameters, (sys.A,sys.B,sys.Bd,sys.C,sys.D)


h = 1/N;

bw = 2;
bphi = 2;
xis = (0:h:1).';

% The simulation uses values \rho=1 and I_\rho=1 (hard-coded in operators
% A, B, C).
rho = 1;
I_rho = 1;

% State variables: 
% x1 = w_z-phi, x2 = rho*w_t
% x3 = phi_z, x4 = I_rho*phi_t

% Boundary conditions:
% x_1(1)=x_3(1)=0
% x_2(0)=x_4(0)=0

% The spatial grids for the approximations of the state variables:
% xi_0 = 0, xi_N = 1
% x1 = x_1(0) ... x_1((N-1)h)
% x3 = x_3(0) ... x_3((N-1)h)
% x2 = x_2(h) ... x_2(Nh)
% x4 = x_4(h) ... x_4(Nh)

ee = ones(N,1);

P012 = 1/h*spdiags([-ee ee],-1:0,N,N);
P021 = 1/h*spdiags([-ee ee],0:1,N,N);

A11 = [sparse(N,N) P012;P021 -bw*speye(N,N)];
A22 = [sparse(N,N) P012;P021 -bphi*speye(N,N)];

P114 = spdiags([-ee 0*ee],-1:0,N,N);
P141 = -P114.';

ZN = sparse(N,N);

A12 = [ZN P114;ZN ZN];
A21 = [ZN ZN;P141 ZN];

A = [A11 A12;A21 A22];


% Control and observation on [4/5 1], both on the state variable x_4(t)
B = [zeros(3*N,1);(xis(2:end)>=.8)];
C = h*B.';


Sys.A = A;
Sys.B = B;
Sys.Bd = sparse(4*N,1);
Sys.C = C;
Sys.D = 0;
Sys.Dd = zeros(size(Sys.C,1),size(Sys.Bd,2));

spgrid = xis;

% Compute the initial state based on w0fun, wd0fun, phi0, and phid0
x0 = [diff(w0fun(xis))/h-phi0fun(xis(1:(end-1)));rho*wd0fun(xis(2:end));diff(phi0fun(xis))/h;I_rho*phid0fun(xis(2:end))];
