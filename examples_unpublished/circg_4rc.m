function [x,y]=circg_4rc(bs,s) 
%CIRCG Geometry File defining the geometry of a circle. 
nrc = 4; %%%%number of rectangles inside the circle
nbs=4*(nrc+1); %define a rectangular inside; 
%nbs=4; %not define a rectangular inside;
rr=1;%rayon
c=[0 0];%center
% set the vertices of the quadrilater
C1=[1/4 1/6]; C2=[-0.6 -1/3]; C3=[-0.4 0.4]; C4 = [0.7 0];
CR=[C1 C2 C3 C4];
rc=[.1 .1 .1 .15 .1 .12 .1 .08];

if nargin==0  
  x=nbs;   
  return 
end
dic=[  0      pi/2   pi       3*pi/2
      pi/2   pi     3*pi/2   2*pi
      1      1      1        1
      0      0      0        0];

dci1=[
      0       0   0       0     
      1       1   1       1
      1       1   1       1  
      2       2   2       2 
    ];

dci2=[
      0       0   0       0     
      1       1   1       1
      1       1   1       1  
      3       3   3       3 
    ];

dci3=[
      0       0   0       0     
      1       1   1       1
      1       1   1       1  
      4       4   4       4 
    ];
dci4=[
      0       0   0       0     
      1       1   1       1
      1       1   1       1  
      5       5   5       5 
    ];
dl=[dic dci1 dci2 dci3 dci4];

if nargin==1   
  x=dl(:,bs);   
  return 
end 

x=zeros(size(s)); 
y=zeros(size(s)); 
[m,n]=size(bs); 
if m==1 & n==1,   
  bs=bs*ones(size(s)); % expand bs 
elseif m~=size(s,1) | n~=size(s,2),   
  error('bs must be scalar or of same size as s'); 
end

if ~isempty(s)
    
% boundary segment 1
ii=find(bs==1);
if length(ii)
x(ii)=rr*cos(s(ii))+c(1);
y(ii)=rr*sin(s(ii))+c(2);
end

% boundary segment 2
ii=find(bs==2);
if length(ii)
x(ii)=rr*cos(s(ii))+c(1);
y(ii)=rr*sin(s(ii))+c(2);
end

% boundary segment 3
ii=find(bs==3);
if length(ii)
x(ii)=rr*cos(s(ii))+c(1);
y(ii)=rr*sin(s(ii))+c(2);
end

% boundary segment 4
ii=find(bs==4);
if length(ii)
x(ii)=rr*cos(s(ii))+c(1);
y(ii)=rr*sin(s(ii))+c(2);
end

for j=0:nrc
    jj=4+4*j;
ii=find(bs==jj+1);
if length(ii)
x(ii)=interp1([dl(1,jj+1),dl(2,jj+1)],[CR(2*j+1)-rc(2*j+1) CR(2*j+1)+rc(2*j+1)],s(ii));
y(ii)=interp1([dl(1,jj+1),dl(2,jj+1)],[CR(2*j+2)-rc(2*j+2) CR(2*j+2)-rc(2*j+2)],s(ii));
end

ii=find(bs==jj+2);
if length(ii)
x(ii)=interp1([dl(1,jj+2),dl(2,jj+2)],[CR(2*j+1)+rc(2*j+1) CR(2*j+1)+rc(2*j+1)],s(ii));
y(ii)=interp1([dl(1,jj+2),dl(2,jj+2)],[CR(2*j+2)-rc(2*j+2) CR(2*j+2)+rc(2*j+2)],s(ii));
end

ii=find(bs==jj+3);
if length(ii)
x(ii)=interp1([dl(1,jj+3),dl(2,jj+3)],[CR(2*j+1)+rc(2*j+1) CR(2*j+1)-rc(2*j+1)],s(ii));
y(ii)=interp1([dl(1,jj+3),dl(2,jj+3)],[CR(2*j+2)+rc(2*j+2) CR(2*j+2)+rc(2*j+2)],s(ii));
end

ii=find(bs==jj+4);
if length(ii)
x(ii)=interp1([dl(1,jj+4),dl(2,jj+4)],[CR(2*j+1)-rc(2*j+1) CR(2*j+1)-rc(2*j+1)],s(ii));
y(ii)=interp1([dl(1,jj+4),dl(2,jj+4)],[CR(2*j+2)+rc(2*j+2) CR(2*j+2)-rc(2*j+2)],s(ii));
end   

end
end
end