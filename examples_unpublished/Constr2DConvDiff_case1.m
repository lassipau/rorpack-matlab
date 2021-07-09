function [x0,Sys] = Constr2DConvDiff_case1(h_param,refine_times,testcase,x0fun)
% Construct the FEM approximation for the 2D reaction-convection-diffusion
% equation "Case 1". Taken from the article "Reduced order controller
% design for robust output regulation" by Lassi Paunonen and Duy Phan, IEEE
% TAC 2020.
%
% h_param = initial mesh size parameter
% refine_times = a number of refinements of the mesh around the actuation
%                and measurement areas
% testcase = determines the "case" of the physical parameters
% x0fun = A function determing the initial state of the system. x0fun is a
%         function of two variables (x0fun = @(x1,x2) ...) and accepts
%         vector inputs
%
% Written by Duy Phan, 2020






% Initiate mesh on the circle, refine near input and output locations
[p, e, t] = refinecircle(h_param,refine_times);
pdeplot(p,e,t);
drawnow


% Construct the FEM approximation fundamental matrices
[fMa,fSt,fGx,fGy]=fullEssentialMatrices(p,e,t);
[St_ii,Ma_ii,Gx_ii,Gy_ii,p_i,Per] = ReorderEssentialMatrices(p,e,fMa,fSt,fGx,fGy);
ni = length(p_i);
nb = length(p) - ni;


% Set the physical parameters of the model
[fa,fb,fc1,fc2] = KnownFunction(testcase);

% Compute the matrices A, B, C of the Galerkin approximation
Sys.A =  ApproxOptA(Ma_ii,St_ii,Gx_ii,Gy_ii, fa,fb,fc1, fc2,p_i);
Sys.B = ApproxOptB(Ma_ii,p_i);
Sys.C = ApproxOptC(Ma_ii,p_i);
Sys.D = zeros(2,2);


x0 = x0fun(p_i(1,:),p_i(2,:)).';



end



function [A] =  ApproxOptA(M,S,Gx,Gy, fa,fb,fc1, fc2, p)
   px = p(1,:); 
   py = p(2,:);
   n = length(px);
   Da = diag(fa(px,py));
   Db = diag(fb(px,py));
   Dc1 = diag(fc1(px,py));
   Dc2 = diag(fc2(px,py)); 
   invM = M\eye(n); 
   A = -S*Da + (M*Db)/2 + (Db*M)/2 + Gx*Dc1 + Gy*Dc2; 
   A = invM*A; 
end

function [B] = ApproxOptB(M,p)
    px = p(1,:); py = p(2,:); n = length(p); 
    %%%The place define controls and observations
    C1=[1/4 1/6]; C2=[-0.6 -1/3]; C3=[-0.4 0.4]; C4 = [0.7 0];
    CR=[C1 C2 C3 C4];
    rc=[.1 .1 .1 .15 .1 .12 .1 .08];
    
    %%Rectangle to put the controls and observations
    R_B = [1 4];
    B = zeros(n,length(R_B));  
    for i = 1:2
        place = R_B(i); 
        x_lo = CR(2*place-1) - rc(2*place-1);
        x_up = CR(2*place-1) + rc(2*place-1);
        y_lo = CR(2*place) - rc(2*place);
        y_up = CR(2*place) + rc(2*place);
        B(:,i) = transpose((px > x_lo).* (px < x_up).*(py > y_lo).*(py < y_up)); 
    end 
    %B = M*B;
%       nb = length(p) - length(p_ii); 
%      B = [B; zeros(nb,2)];
%      vv =  Per'*(B(:,1)+B(:,2));
%      pdeplot(p,e,t,'xydata',vv,'zdata',vv,'colormap','jet','mesh','on');

end 



function   [C] = ApproxOptC(M,p)
    px = p(1,:); py = p(2,:); n = length(p); 
    %%%The place define controls and observations
    C1=[1/4 1/6]; C2=[-0.6 -1/3]; C3=[-0.4 0.4]; C4 = [0.7 0];
    CR=[C1 C2 C3 C4];
    rc=[.1 .1 .1 .15 .1 .12 .1 .08];
    
    %%Rectangle to put the controls and observations
    R_C = [2 3];
    C = zeros(length(R_C),n);  
    for i = 1:length(R_C)
        place = R_C(i); 
        x_lo = CR(2*place-1) - rc(2*place-1);
        x_up = CR(2*place-1) + rc(2*place-1);
        y_lo = CR(2*place) - rc(2*place);
        y_up = CR(2*place) + rc(2*place);
        C(i,:) = ((px > x_lo).* (px < x_up).*(py > y_lo).*(py < y_up)); 
    end 
    C = C*M;
end


function [fa,fb,fc1,fc2] = KnownFunction(testcase)

    switch testcase
        case 1 
            fa = @(x,y) 0.5*(x.^2 + 0.1);
            fb = @(x,y) 15*x;
            i = 1; j = 2; k = 3; l = 4;
            fc1 =  @(x,y) cos(i*x) - sin(j*y);
            fc2 = @(x,y) sin(k*x) + cos(l*y);
        case 2 
             fa = @(x,y) (x.^2 + 1);
            fb = @(x,y) 15;
            i = 1; j = 2; k = 3; l = 4;
            fc1 =  @(x,y) cos(i*x) - sin(j*y)+5;
            fc2 = @(x,y) sin(k*x) + cos(l*y)-5;
        case 3
            fa = @(x,y) 0.25;
            fb = @(x,y) 10*(x+y);
             i = 1; j = 2; k = 3; l = 4;
            fc1 =  @(x,y) cos(i*x) - sin(j*y);
            fc2 = @(x,y) sin(k*x) + cos(l*y);
         case 4
            fa = @(x,y) 0.5;
            fb = @(x,y) 10;
             i = 1; j = 2; k = 3; l = 4;
            fc1 =  @(x,y) cos(i*x) - sin(j*y);
            fc2 = @(x,y) sin(k*x) + cos(l*y);
%             % LASSI ADDED rescaling
%             fc1 = @(x,y) 0.5*(cos(i*x) - sin(j*y));
%             fc2 = @(x,y) 0.3*(sin(k*x) + cos(l*y));
         case 5 
            fa = @(x,y) (x.^2 + 0.5);
            fb = @(x,y) 10;
            fc1 =  @(x,y) 0;
            fc2 = @(x,y) 0;   
         case 6 
            fa = @(x,y) (x.^2 + 0.5);
            fb = @(x,y) 10*(x+y);
            fc1 =  @(x,y) 0;
            fc2 = @(x,y) 0;   
            
%         case 2 
%         st = 1;     
%         nu = 2;
%         fa = 60+3*sin(xx).^2+5*cos(yy);
%         vb1 = zeros(1,np);
%         vb2 = zeros(1,np);
%         case 3
%         st = 1;
%         nu = 2;
%         fa = 60+3*sin(xx).^2+5*cos(yy);
%         vb1 = 10*sin(xx);
%         vb2 = cos(yy);
%         case 4
%         nu = 3;
%         fa = 60+3*sin(xx).^2+5*cos(yy);
%         vb1 = exp(tt/2)*sin(xx);
%         vb2 = cos(yy);
%         case 5
%         nu = 0.25;
%         i = 1; j = 2; k = 3; l = 4; m = 5; n = 6;
%         fa = sin(tt)*cos(i*xx) - sin(5*tt)*sin(j*yy)+3;
%         vb1 = cos(tt)*sin(k*xx) + cos(3*tt)*cos(l*yy);
%         vb2 = sin(tt)*sin(m*xx) + cos(2*tt)*sin(n*yy);   
%         case 6
%         nu = 0.25;
%         i = 1; j = 1; k = 1; l = 1; m = 1; n = 1;
%         fa =  sin(tt)*cos(i*xx) - sin(5*tt)*sin(j*yy)+3;
%         vb1 = cos(tt)*sin(k*xx) + cos(3*tt)*cos(l*yy);
%         vb2 = sin(tt)*sin(m*xx) + cos(2*tt)*sin(n*yy);
%         case 7
%         nu = 0.1;
%         nu = 0.25;
%         i = 1; j = 2; k = 2; l = 1; m = 1; n = 1;
%         fa =  (sin(tt)*cos(i*xx) - sin(5*tt)*sin(j*yy))+3;
%         vb1 = cos(tt)*sin(k*xx) + cos(3*tt)*cos(l*yy);
%         vb2 = sin(tt)*sin(m*xx) + cos(2*tt)*sin(n*yy);
%         case 8
%         nu = 0.25;
%         i = -1; j = 5; k = 3; l = 1; m = 1; n = 5;
%         fa =  (sin(tt)*cos(i*xx) - sin(5*tt)*sin(j*yy))+3;
%         vb1 = cos(tt)*sin(k*xx) + cos(3*tt)*cos(l*yy);
%         vb2 = sin(tt)*sin(m*xx) + cos(2*tt)*sin(n*yy);
%         case 9
%         st = 1;
%         nu = 0.5; 
%         fa = 10*ones(1,np); 
%         vb1 = zeros(1,np);
%         vb2 = zeros(1,np); 
%         case 10
%         st = 1; 
%         nu = 0.5; 
%         fa = 10*ones(1,np); 
%         vb1 = xx.^2;
%         vb2 = sin(yy); 
%         case 11
%         nu = 0.1;
%         i = 6; j = -2; k = 5; l = 3; m = 4; n = 1;
%         fa =  sin(tt)*cos(i*xx) - sin(5*tt)*sin(j*yy)+2;
%         vb1 = cos(tt)*sin(k*xx) + cos(3*tt)*cos(l*yy);
%         vb2 = sin(tt)*sin(m*xx) + cos(2*tt)*sin(n*yy); 
%         case 12
%         st = 1; 
%         nu = 1; 
%         fa = 10*ones(1,np)-2*xx-cos(yy); 
%         vb1 = xx.^2;
%         vb2 = sin(yy);
%         case 13
%         st = 1; 
%         nu = 0.5; 
%         fa = 10*ones(1,np)-2*xx-cos(yy); 
%         vb1 = xx.^2;
%         vb2 = sin(yy);
%         case 14 
%         nu = 0.25; 
%         fa = 10*ones(1,np)-2*xx-sin(2*tt)*cos(yy); 
%         vb1 = cos(3*tt)*xx.^3;
%         vb2 = sin(tt)*sin(yy);
    end
    
end
   

function [fMa,fSt,fDx,fDy]=fullEssentialMatrices(p,e,t)
%
% ---------- In the plane ------------
% Given the mesh (p,e,t); constructs the stiffness and rhs matrix
% and the boundary mass matrix
nzmax=size(p,2)^2;
npt=size(p,2);
ned=size(e,2);
ntr=size(t,2);

fprintf('\n npt: %2.0f\n ned: %2.0f\n ntr: %2.0f\n\n', npt,ned,ntr)
fprintf('Stiffness and Mass:  ')
tic
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     STIFFNESS and MASS matrices
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% for each triangle we find the respective contributions
%
%


fMa=sparse(npt,npt);
fSt=sparse(npt,npt);
fDx=sparse(npt,npt);
fDy=sparse(npt,npt);
%maxtp=zeros(1,ne);
% indicator that point e(1,l) is a vertice of t(k)
pt=zeros(ned,ntr);
for k=1:ntr
    vk=t(1:3,k);
    v12=p(:,vk(2))-p(:,vk(1));
    v13=p(:,vk(3))-p(:,vk(1));
    % Let fi be the linear test function with value 1
    % at vertice i and value 0 at the others:
    % -------------STIFFNESS---------------
    %     df1(v12)=-1=df1(v13)
    %     df2(v12)=1=df3(v13)
    %     df2(v13)=0=df3(v12)
    % Computing [df1/dx; df1/dy] we denote by df(i,:)
    %A=[v12 v13];
    areatk=(v12(1)*v13(2)-v12(2)*v13(1))/2;
    % detA=2*areatk
    invA=1/(2*areatk)*[ v13(2) -v13(1) ; -v12(2) v12(1) ];
    D1=invA*[1 ; 0]; D2=invA*[0 ; 1];
    df(1,:)=[-D1(1)-D1(2) -D2(1)-D2(2)];
    df(2,:)=[D1(1) D2(1)];
    df(3,:)=[D1(2) D2(2)];
    sc=zeros(3,3);
    for i=1:3
        for j=i:3
            int=df(i,:).*df(j,:);
            sc(i,j)=int(1)+int(2);
        end
    end
    ii=[vk(1) vk(1) vk(1) vk(2) vk(2) vk(2) vk(3) vk(3) vk(3)]';
    jj=[vk(1) vk(2) vk(3) vk(1) vk(2) vk(3) vk(1) vk(2) vk(3)]';
    ss=areatk*[sc(1,1) sc(1,2) sc(1,3) sc(1,2) sc(2,2) sc(2,3) sc(1,3) sc(2,3) sc(3,3)]';
    lkstM=sparse(ii,jj,ss,npt,npt,nzmax);%localk
    fSt=fSt+lkstM;
    %  -------------MASS------------------
    % We can check that
    %         f1(a*v12+b*v13)=1-a-b
    %         f2(a*v12+b*v13)=a
    %         f3(a*v12+b*v13)=b
    % and, the area element in this triangle is
    %            |v12||v13|sin(angle(v12,v13)).da^db=2*areatk.da^db
    % easily we find that the ... given by any triangle tk to the rhs
    % matrix at nodes (ti,tj) is mM(i,i)=1/6*areatk, and
    % mM(i,j)=1/12*areatk if i~=j.
    sm=1/12*areatk*[2 1 1 1 2 1 1 1 2]';
    lkrhM=sparse(ii,jj,sm,npt,npt,nzmax);%localk
    fMa=fMa+lkrhM;
    %--------------Dx&Dy matrices------------
    ssDx=1/6*[(p(2,vk(2))-p(2,vk(3)))*ones(1,3); (p(2,vk(3))-p(2,vk(1)))*ones(1,3); (p(2,vk(1))-p(2,vk(2)))*ones(1,3) ]';
    sx =[ssDx(1,:) ssDx(2,:) ssDx(3,:)];
    lkDx=sparse(ii,jj,sx,npt,npt,nzmax);%localk
    fDx=fDx+lkDx;
    ssDy=1/6*[(p(1,vk(3))-p(1,vk(2)))*ones(1,3); (p(1,vk(1))-p(1,vk(3)))*ones(1,3); (p(1,vk(2))-p(1,vk(1)))*ones(1,3)]';
    sy =[ssDy(1,:) ssDy(2,:) ssDy(3,:)];
    lkDy=sparse(ii,jj,sy,npt,npt,nzmax);%localk
    fDy=fDy+lkDy;
    
    %----------
    %     % The following 'for' is needed only if not Dirichlet b.c.
    %     for l=1:ned
    %         %check if t intersects the origin of some edge (on the boundary)
    %         c1=t(1,k)==e(1,l)|t(2,k)==e(1,l)|t(3,k)==e(1,l);
    %         %c2=t(1,k)==e(2,l)|t(2,k)==e(2,l)|t(3,k)==e(2,l);
    %         if c1==1
    %             %ptb(l)=e(1,l);
    %             %maxtp(l)=maxtp(l)+1;
    %             pt(l,k)=1;%indicator that point e(1,l) is a vertice of tk
    %         end
    %     end
    %----------so, later maybe try to avoid its comp. if no need)
end

toc
end

function [St_ii,Ma_ii,Dx_ii,Dy_ii,p_ii,Per] = ReorderEssentialMatrices(p,e,fMa,fSt,fDx,fDy)
npt=length(p);
ned=length(e);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%        Boundary points last

bpts=sparse(1,npt);
ipts=1:1:npt;
nb=0; %% number of boundary points
% we have to know our geometry in this case
% e[6 7] in {[1 0],[0 3]} correspond to
% the real boundary
for i=1:ned
    vi=find(ipts==e(1,i));
    if (isequal(e(6:7,i),[1 0]') || isequal(e(6:7,i),[0 3]')) % [1 0] -- exterior boundary,  [0 1] -- for hole(s)
        bpts(vi)=1;
        nb=nb+1;
    end
end

[opts,Iib]=sortrows(bpts');


% % % permutation matrix
nmax=npt^2;
ii=linspace(1,npt,npt)';
jj=Iib;
vv=ones(npt,1);
Per=sparse(ii,jj,vv,npt,npt,nmax);
px = Per*p(1,:)';
py = Per*p(2,:)';
p_ii = [px(1:npt-nb), py(1:npt-nb)];
p_ii = transpose(p_ii);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%          FINAL MATRICES
%

% full matrices in new coordinates (i.e., boundary points last)
foSt=Per*fSt*Per';
foMa=Per*fMa*Per';
foDx=Per*fDx*Per';
foDy=Per*fDy*Per';

%%%% Mass in Blocks

foMa_ii=[eye(npt-nb)  zeros(npt-nb,nb)]*foMa*[eye(npt-nb); zeros(nb,npt-nb)];
%foMa_bi=[zeros(nb,npt-nb) eye(nb)]*foMa*[eye(npt-nb); zeros(nb,npt-nb)];

%foMa_ib=[eye(npt-nb) zeros(npt-nb,nb)]*foMa*[zeros(npt-nb,nb) ; eye(nb)];
%foMa_bb=[zeros(nb,npt-nb)  eye(nb)]*foMa*[zeros(npt-nb,nb) ; eye(nb)];

%%%% Stiffness in Blocks


foSt_ii=[eye(npt-nb)  zeros(npt-nb,nb)]*foSt*[eye(npt-nb) ; zeros(nb,npt-nb)];
%foSt_bi=[zeros(nb,npt-nb) eye(nb)]*foSt*[eye(npt-nb) ; zeros(nb,npt-nb)];

%foSt_ib=[eye(npt-nb)  zeros(npt-nb,nb)]*foSt*[zeros(npt-nb,nb) ; eye(nb)];
%foSt_bb=[zeros(nb,npt-nb)  eye(nb)]*foSt*[zeros(npt-nb,nb) ; eye(nb)];



%%%% Stiffness in Blocks

foSt_ii=[eye(npt-nb)  zeros(npt-nb,nb)]*foSt*[eye(npt-nb) ; zeros(nb,npt-nb)];
%foSt_bi=[zeros(nb,npt-nb) eye(nb)]*foSt*[eye(npt-nb) ; zeros(nb,npt-nb)];

%foSt_ib=[eye(npt-nb)  zeros(npt-nb,nb)]*foSt*[zeros(npt-nb,nb) ; eye(nb)];
%foSt_bb=[zeros(nb,npt-nb)  eye(nb)]*foSt*[zeros(npt-nb,nb) ; eye(nb)];


%%%% Dx, Dy in Blocks

foDx_ii=[eye(npt-nb)  zeros(npt-nb,nb)]*foDx*[eye(npt-nb); zeros(nb,npt-nb)];
%foDx_bi=[zeros(nb,npt-nb) eye(nb)]*foDx*[eye(npt-nb); zeros(nb,npt-nb)];

%foDx_ib=[eye(npt-nb) zeros(npt-nb,nb)]*foDx*[zeros(npt-nb,nb) ; eye(nb)];
%foDx_bb=[zeros(nb,npt-nb)  eye(nb)]*foDx*[zeros(npt-nb,nb) ; eye(nb)];

foDy_ii=[eye(npt-nb)  zeros(npt-nb,nb)]*foDy*[eye(npt-nb); zeros(nb,npt-nb)];
%foDy_bi=[zeros(nb,npt-nb) eye(nb)]*foDx*[eye(npt-nb); zeros(nb,npt-nb)];

%foDy_ib=[eye(npt-nb) zeros(npt-nb,nb)]*foDy*[zeros(npt-nb,nb) ; eye(nb)];
%foDy_bb=[zeros(nb,npt-nb)  eye(nb)]*foDx*[zeros(npt-nb,nb) ; eye(nb)];

% %%%%Control Q in blocks
% for j=1:nt
% Q_i(:,:,j)= [eye(npt-nb)  zeros(npt-nb,nb)]*Q(:,:,j)*[eye(npt-nb); zeros(nb,npt-nb)];
% Q_b(:,:,j)=[eye(npt-nb) zeros(npt-nb,nb)]*Q(:,:,j)*[zeros(npt-nb,nb) ; eye(nb)];
% end

St_ii=foSt_ii;
% St_ib=foSt_ib;

Ma_ii=foMa_ii;
% Ma_ib=foMa_ib;

Dx_ii=foDx_ii;
% Dx_ib=foDx_ib;

Dy_ii=foDy_ii;
% Dy_ib=foDy_ib;
end


function [p, e, t] = refinecircle(h,ii)
%%%%Updata 25 June 2018
%%% The role of function is to create a mesh where the maximum
%%% length of triangle's sides is h. Then, refine 3 sub-rectangles.
%%%ii: times refining mesh inside the rectangle.
%%%oo: times refining mesh outside the rectangle.
fprintf('\t Create the mesh and Refine it in the sub-rectangle(s) \n')
[p, e, t]= initmesh('circg_4rc','hmax',h);
for jj = 1:ii
    ptemp = p;
    etemp = e;
    ttemp = t;
    p = []; e =[]; t=[]; tri = [];
    for ind = 1:length(ttemp)
        if ttemp(4,ind) ~= 1
            tri = [tri; ind];
        end
    end
    [p,e,t] = refinemesh('circg_4rc',ptemp,etemp,ttemp,tri,'regular');
    ptemp = []; etemp =[]; ttemp=[];
end
%     for jj = 1:oo
%         ptemp = p;
%         etemp = e;
%         ttemp = t;
%         p = []; e =[]; t=[];
%         [tri] = FindTriInRect(ptemp,etemp,ttemp,0-h/2,1/2+h/2,0-h/2,1/3+h/2);
%         triout = setdiff(1:length(ttemp),tri);
%         [p,e,t] = refinemesh('circg',ptemp,etemp,ttemp,triout','regular');
%         ptemp = []; etemp =[]; ttemp=[];
%     end
end

