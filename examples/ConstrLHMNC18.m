function [x0,spgrid,sys] = ConstrLHMNC18(x0fun,N)
% [x0,spgrid,sys] = ConstrHeat2Dtest1(c,x0fun,N)
% Usage: Can also use solely for defining the initial state x0
% x0 = initial state
% spgrid = spatial grid of the domain
% sys = system parameters, (sys.A,sys.B,sys.Bd,sys.C,sys.D)


h = 1/N;

bw = 2;
bphi = 2;
xis = (0:h:1).';


% xi_0 = 0, xi_N = 1
% x1 = x_1(0) ... x_1((N-1)h)
% x3 = x_3(0) ... x_3((N-1)h)
% x2 = x_2(h) ... x_2(Nh)
% x4 = x_4(h) ... x_4(Nh)

% x_1(1)=x_3(1)=0
% x_2(0)=x_4(0)=0

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


% Control and observation on [4/5 1], both on the element x4

B = [zeros(3*N,1);(xis(2:end)>=.8)];
C = h*B.';


sys.A = A;
sys.B = B;
sys.Bd = sparse(4*N,1);
sys.C = C;
sys.D = 0;

spgrid = xis;
x0 = x0fun(xis);

