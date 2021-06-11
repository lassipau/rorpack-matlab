function [Dop,spgrid] = DiffOp1d(cfun,spgrid,BCtype)

% Sparse discrete 1D diffusion operator A defined by (Au)(x)=(c(x)*u'(x))' 
% with Dirichlet or Neumann boundary conditions.
% 
%     Parameters:
%     ---------
%     spgrid : 
%         The discretized (uniform) spatial grid on the interval [0,L] so that
%         spgrid(1)=0 and spgrid(end)=L, 
%     cfun : function
%         Defines the spatially dependent coefficient function c(x)
%     BCtype : string
%         Boundary condition configuration, either
%         'DD' for Dirichlet at x=0, Dirichlet at x=1
%         'DN' for Dirichlet at x=0, Neumann at x=1
%         'ND' for Neumann at x=0, Dirichlet at x=1
%         'NN' for Neumann at x=0, Neumann at x=1
% 
%     Returns
%     ---------
%     Dop :
%         The square matrix representing the diffusion operator.
%     spgrid :
%         An adjusted spatial grid:

N = size(spgrid,2);
assert(N > 0);
h = 1/(N-1);

% We use the approximation:  
% (c(x)u'(x))'(x) ~ 1/h^2*(c(x+h/2)*(u(x+h)-u(x))-c(x-h/2)*(u(x)-u(x-h))

% points x_k+h/2
spmidpoints = 0.5*(spgrid(1:N-1)+spgrid(2:N));
cmid = cfun(spmidpoints);

% Diffusion operator base (Neumann-Neumann case)
cminus = [cmid 0].';
cplus = [0 cmid].';
Dop = spdiags([cminus -(cminus+cplus) cplus],-1:1,N,N);
Dop(1,1) = -(cfun(spgrid(1))+cfun(spgrid(2)));
Dop(1,2) = cfun(spgrid(1))+cfun(spgrid(2));
Dop(N,N) = -(cfun(spgrid(N))+cfun(spgrid(N-1)));
Dop(N,N-1) = cfun(spgrid(N))+cfun(spgrid(N-1));
Dop = 1/h^2*Dop;

% Neumann-Neumann BCs
if isequal(BCtype,'NN')
    % No need to change the grid
    spgrid = spgrid;
elseif isequal(BCtype,'ND')
    spgrid = spgrid(1:N-1);
    Dop = Dop(1:N-1,1:N-1);
elseif isequal(BCtype,'DN')
    spgrid = spgrid(2:N);
    Dop = Dop(2:N,2:N);
elseif isequal(BCtype,'DD')
    spgrid = spgrid(2:N-1);
    Dop = Dop(2:N-1,2:N-1);
else
    error('Unrecognised boundary condition types.')
end