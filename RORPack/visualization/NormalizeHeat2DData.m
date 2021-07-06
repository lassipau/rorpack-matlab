function normalized2DData = NormalizeHeat2DData(N, data, xgrid, ygrid)
% Normalize the data (solutions of the controlled PDE) used in plotting and
% animating the state of the 2D heat equation

nn = 0:N-1;
mm = 0:N-1;
te = size(data,2);
result = zeros(size(xgrid,1),size(xgrid,2),te);

for i = 1:te
  gat = reshape(data(1:N*N, i), [N N]);
  for j = 0:N*N-1
    k1 = floor(j / N)+1;
    k2 = rem(j,N)+1;
    result(:, :, i) = result(:, :, i) + gat(k1, k2) .* phinm(xgrid, ygrid, nn(k1), mm(k2));
  end
end

normalized2DData = reshape(result, [size(xgrid,1)^2 te]);

% Computes eigenfunctions of A.
% phi_00=1
% phi_n0=sqrt(2)*cos(pi*n*x1) for n>=1
% phi_0m=sqrt(2)*cos(pi*m*x2) for m>=1
% phi_nm=sqrt(2)*cos(pi*n*x1)*cos(pi*m*x2) for n,m>=1
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