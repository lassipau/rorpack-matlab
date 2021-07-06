function [x0,spgrid,Sys] = ConstrHeat2DCase3(cval,x0fun,N)
% Construct the system operators for the 2D heat equation example "Case 3".
% The construction is the same as for the "Case 1", but the system is not
% stabilized with output feedback.

xx = linspace(0,1,N);
yy = xx;
[xxg,yyg] = meshgrid(xx,yy);
spgrid.xx = xxg;
spgrid.yy = yyg;

x0 = reshape(x0fun(xxg,yyg),N^2,1);

% Construct the system using modal approximation
M = N;
nn = [0:N-1].';
mm = [0:M-1];
nnmm = nn*ones(1,M) + ones(N,1)*mm;
glnm = reshape(-nnmm.^2 * pi^2,1,N*M)';
A0 = diag(glnm);

b11 = 0.5;
b12 = (1/sqrt(2)) * ones(1,M-1);
b13 = sqrt(2)*sin(nn(2:size(nn,1),:)*(pi/2))./(nn(2:size(nn,1),:)*pi);
b14 = (2*sin(nn(2:size(nn,1),:)*(pi/2))./(nn(2:size(nn,1),:)*pi)) * (ones(1,M-1));
b1 = [b11 b12;b13 b14];

b21 = 0.5;
b22 = (1/sqrt(2)) * (-1) * ones(1,M-1);
b23 = -sqrt(2)*sin(nn(2:size(nn,1),:)*(pi/2))./(nn(2:size(nn,1),:)*pi);
b24part = ones(1,M-1);
b24part(2:2:end) = -b24part(2:2:end);
b24 = (2*sin(nn(2:size(nn,1),:)*(pi/2))./(nn(2:size(nn,1),:)*pi)) * b24part;
b2 = [b21 b22;b23 b24];

Sys.A = A0;
B = [b1(:) b2(:)];
Sys.B = B;
C = 2*B';
Sys.C = C;
Sys.D = zeros(size(C,1),size(B,2));
Sys.Bd = zeros(N*M,1);
Sys.Dd = zeros(size(C,1),1);

end