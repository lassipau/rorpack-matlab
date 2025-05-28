function [x0,spgrid,Sys] = ConstrHeat2DCase1(cval,x0fun,N)
% Construct the system operators for the 2D heat equation example "Case 1".
% The system has 2 boundary inputs and 2 collocated boundary outputs, and
% an additional boundary disturbance input. The system is impedance
% passive, and it is prestabilized with negative output feedback. The
% numerical approximation is based on the modal approximation of the 2D
% Laplacian.

xx = linspace(0,1,N);
yy = xx;
[xxg,yyg] = meshgrid(xx,yy);
spgrid.xx = xxg;
spgrid.yy = yyg;


% Construct the system using modal approximation
M = N;
nn = (0:N-1).';
mm = (0:M-1);
nnmm = nn*ones(1,M) + ones(N,1)*mm;

% x0 = reshape(x0fun(xxg,yyg),N^2,1);
x0 = zeros(N);
for ind1 = 1:N
    for ind2 = 1:M
        x0(ind1,ind2) = integral2(@(x,y) phinm(x,y,ind1-1,ind2-1).*x0fun(x,y) ,0,1,0,1);
    end
end
x0 = x0(:);

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

B = [b1(:) b2(:)];
Sys.B = B;
C = B';
Sys.C = C;
Sys.D = zeros(size(C,1),size(B,2));
Sys.Bd = zeros(N*M,1);
Sys.Dd = zeros(size(C,1),1);

% Exponential prestabilization with negative output feedback.
Sys.A = A0 - B*C;

end

function eigFun = phinm(x1, x2, n, m)
if n == 0 && m == 0
    eigFun = ones(size(x1));
end
if n == 0 && m >= 1
    eigFun = sqrt(2) * cos(pi * m * x2);
end
if n >= 1 && m == 0
    eigFun = sqrt(2) * cos(pi * n * x1);
end
if n >= 1 && m >= 1
    eigFun = 2 * cos(pi * n * x1) .* cos(pi * m * x2);
end
end