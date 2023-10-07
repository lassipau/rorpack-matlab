function [x0,Sys,phin] = ConstrWave1DCase3(dfun,xi0,w0fun,wd0fun,N) 
% Modal approximation of a 1D Wave equation with viscous damping profile
% determined by 'dfun'.
% Neumann boundary input (at xi=0) and Dirichlet boundary condition w(1,t)=0.
% Measured outputs at a specified point xi0, i.e.,
% y(t)=(w(xi0,t),w_t(xi0,t))^T. Input disturbance at xi=0.
%
% Usage: Can also be used solely for defining the initial state x0
% w0fun = initial displacement, function handle
% wd0fun = initial velocity, function handle
% xi0 = the x-coordinate of the location of the measurements
% sys = system parameters, (sys.A,sys.B,sys.Bd,sys.C,sys.D)
% phin = Normalized eigenfunction of the Laplacian


phin = @(xi,n) sqrt(2)*cos((n+1/2)*pi*xi);
phidn = @(xi,n) -(n+1/2).*pi.*sqrt(2).*sin((n+1/2)*pi*xi);

% The system operator A is constructed has the form
% A=[0,I;Lap,-D], where 'Lap' is the approximate Laplacian,
% D>=0 is the damping operator,
% The input and output operators have the forms
% B = [0;B0], C = [C1,0;0,C2].
% The initial state has the form x0 = [x01;x02]

nn = 0:(N-1);
Lap = diag(-(nn+1/2).^2*pi^2);
B0 = phin(0,nn');
C1 = phin(xi0,nn);
C2 = phin(xi0,nn);


D = zeros(N);
x01 = zeros(N,1);
x02 = zeros(N,1);
for indn = 0:(N-1)
    for indm = 0:(N-1)
        D(indm+1,indn+1) = integral(@(xi) dfun(xi).*phin(xi,indn).*conj(phin(xi,indm)),0,1);
    end
    x01(indn+1) = integral(@(xi) w0fun(xi).*conj(phin(xi,indn)),0,1);
    x02(indn+1) = integral(@(xi) wd0fun(xi).*conj(phin(xi,indn)),0,1);
end


A = [zeros(N),eye(N);Lap,-D];
B = [zeros(N,1);B0];
C = [C1,zeros(1,N);zeros(1,N),C2];

x0 = [x01;x02];


Sys.A = sparse(A);
Sys.B = B;
Sys.Bd = B; % Input disturbance
Sys.C = C;
Sys.D = [0;0];
Sys.Dd = [0;0];



