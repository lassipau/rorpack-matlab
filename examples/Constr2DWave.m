function [Sys,svals,sfvals,phinm1,phinm2,phiRnm,psiTHm1,psiTHm2,Q] = Constr2DWave(Nvals,Mvals)

% eigenvalues and -functions
% m = index of the angular mode
% n = index of the radial mode (dep on m)
zerfun = @(s,m) besselj(m,s).*(m/2.*bessely(m,2*s)-s.*bessely(m+1,2*s)) ...
                -(m/2*besselj(m,2*s)-s.*besselj(m+1,2*s)).*bessely(m,s);

svals = zeros(Mvals,Nvals);
sfvals = zeros(Mvals,Nvals);

for indm = 0:(Mvals-1)
  if indm <= 2
    frstzer = 3;
  elseif indm <=5
    frstzer = 4;
  elseif indm <= 9
    frstzer = 5.5;
  elseif indm <= 11
    frstzer = 6.5;
  elseif indm <= 13
    frstzer = 7.5;
  elseif indm <= 16
    frstzer = 10.5;
  elseif indm <= 21
    frstzer = 12;
  elseif indm <= 26
    frstzer = 14.5;
  elseif indm <= 28
    frstzer = 15.6;
  else
    frstzer = 18.5;
  end
  [svals(indm+1,1),sfvals(indm+1,1)] = ...
      fzero(@(s) zerfun(s,indm),[1.1 frstzer]);
  for indn = 2:Nvals
    [svals(indm+1,indn),sfvals(indm+1,indn)] = ...
      fzero(@(s) zerfun(s,indm),svals(indm+1,indn-1)+[1,4]);
  end
end

for indm = 1:Mvals
  if length(round(svals(indm,:),2))>length(unique(round(svals(indm,:),2)))
    error('multiple roots in eigenvalue search!')
  end
end

% unnormalised radial  basis funtions
phinbas = @(r,m,snm) besselj(m,snm*r).* ...
    bessely(m,snm)-besselj(m,snm).*bessely(m,snm*r);

phinorms = zeros(Mvals,Nvals);
Jcoeffs = phinorms; 
Ycoeffs = phinorms;

for indm = 0:(Mvals-1)
  for indn = 1:Nvals
    phinorms(indm+1,indn) = integral(@(r) r.*phinbas(r, ...
        indm,svals(indm+1,indn)).^2,1,2,'Abstol',1e-6,'Reltol',1e-12);
    Jcoeffs(indm+1,indn) = besselj(indm,svals(indm+1,indn));
    Ycoeffs(indm+1,indn) = bessely(indm,svals(indm+1,indn));
  end
end

% Radial component only (needed in defining B and C)
phiRnm = @(r,n,m) phinorms(m+1,n)^(-1/2)*phinbas(r,m,svals(m+1,n));
psiTHm1 = @(theta,m) 1/sqrt(pi).*cos(m*theta).*(m>0) ...
    +ones(size(theta))./sqrt(2*pi).*double(m==0);
psiTHm2 = @(theta,m) 1/sqrt(pi).*sin(m*theta).*(m>0);

% Cosine basis for theta, for m>=0
phinm1 = @(r,theta,n,m) phiRnm(r,n,m).*psiTHm1(theta,m);
% Sine basis for theta, only for m>0
phinm2 = @(r,theta,n,m) phiRnm(r,n,m).*psiTHm2(theta,m);

% approximate operators A, B, C
% state vector structure: 2x2 blocks associated to modes
% n=1, m=0..M-1, n=2, m=0..M-1, ... n=N, m=0..M-1

% Dimension of the input/output space = 2*Mvals-1 (Mvals*cosines +
% (Mvals-1)*sines)

A = zeros(2*Nvals*(2*Mvals-1));
B = zeros(2*Nvals*(2*Mvals-1),2*Mvals-1);
C = zeros(2*Mvals-1,2*Nvals*(2*Mvals-1));

% Each m gets has two 2x2 blocks, one for cosine and one for sine, except
% m=0, which only gets one 2x2 block
for indn = 1:Nvals
  for indm = 0:(Mvals-1)
    if indm==0
      indran = 2*(indn-1)*(2*Mvals-1) + (1:2);
      A(indran,indran) = [0 1;-svals(indm+1,indn)^2 0];
      B(indran,1) = [0;2*conj(phiRnm(2,indn,0))];
      C(1,indran) = [0 phiRnm(2,indn,0)];
    else
      indran = 2*(indn-1)*(2*Mvals-1) + 4*indm - 2 + (1:4);
      A(indran,indran) = blkdiag([0 1;-svals(indm+1,indn)^2 0],[0 1; ...
          -svals(indm+1,indn)^2 0]);
      B(indran,2*indm+[0 1]) = blkdiag([0;2*conj(phiRnm(2,indn,indm))], ...
          [0;2*conj(phiRnm(2,indn,indm))]);
      C(2*indm+[0 1],indran) = blkdiag([0 phiRnm(2,indn,indm)], ...
          [0 phiRnm(2,indn,indm)]);
    end
  end
end
Q = 3; % output feedback gain
As = A-Q*B*C; % stabilized operator A

Sys.A = As;
Sys.B = B;
Sys.C = C;