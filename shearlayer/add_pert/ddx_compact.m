function [RHS]=ddx_compact(f,Lx,nx,ny,nz)
dx = Lx/nx;  

alpha06d1=1/3;             a06d1=(14/9)/2;             b06d1=(1/9)/4;  %first der coeff

%alpha=3;      p=17/6;       q=3/2;       r=3/2;        s=-1/6;  %edge nodes

%qhat=a06d1;   rhat=b06d1;  
%alpha_hat=alpha06d1; %interior nodes

%q_p=3/4;      alpha_p=1/4;  %node2

%alpha_pp = ((40*alpha_hat-1)*q+ 7*(4*alpha_hat-1)*s)/(16*(alpha_hat + 2)*q + 8*(1-4*alpha_hat)*s);
%q_pp= (1/3)*(alpha_pp + 2);r_pp= (1/12)*(4*alpha_pp-1);s_pp=0;  %node-3

%w1 = (2*alpha_hat + 1)/(2*(q + s));        
%w2 =((8*alpha_hat + 7)*q - 6*(2*alpha_hat + 1)*r+ (8*alpha_hat + 7)*s)/(9*(q + s));
%w3 = (4*(alpha_hat + 2)*q + 2*(1 - 4*alpha_hat)*s)/(9*(q + s));     %weights

% Periodic for now

a06 = a06d1/dx;
b06 = b06d1/dx;

RHS=zeros(nx,ny,nz);
RHS(1,:,:) =      a06 *(f(2,:,:)- f(nx,:,:))        + b06*(f(3,:,:)-f(nx-1,:,:));
RHS(2,:,:) =      a06 *(f(3,:,:)- f(1,:,:))         + b06*(f(4,:,:)-f(nx,:,:));
RHS(3:nx-2,:,:) = a06*(f(4:nx-1,:,:)-f(2:nx-3,:,:)) + b06*(f(5:nx,:,:)-f(1:nx-4,:,:));
RHS(nx-1,:,:) =   a06*( f(nx,:,:)-f(nx-2,:,:))      + b06*(f(1,:,:)-f(nx-3,:,:));
RHS(nx,:,:)   =   a06 *(f(1,:,:)- f(nx-1,:,:))      + b06*( f(2,:,:)- f(nx-2,:,:) );


LU = zeros(nx,5);
bc=LU(:,1);     h=LU(:,2); c=LU(:,3); aa=LU(:,4); v=LU(:,5);
b=alpha06d1;d=1;a=alpha06d1;

c(1,1) = d;
v(1,1) = b;
h(1,1) = a/c(1,1);

% Step 1
for i = 2:nx-1
c(i,1) = d - ( b/c(i-1,1) ).*a;
end 
for i = 2:nx-2
v(i,1) = -( b/c(i-1,1) ).*v(i-1,1);
h(i,1) = -( a/c(i,1) ).*h(i-1,1);
end 
v(nx-1,1) = a - ( b/c(nx-2,1) ).*v(nx-2,1);
h(nx-1,1) = ( b - h(nx-2,1).*a ) ./ c(nx-1,1);
c(nx,1)   = d - sum( h(1:nx-1,1).*v(1:nx-1,1) );

bc(2:nx-1,1) = b ./ c(1:nx-2,1);
aa(1:nx-2,1) = a;
c(:,1) = 1./c(:,1);

% Overwrite aa by aa*c
aa(:,1) = aa(:,1).*c(:,1);

%Overwrite v by v*c
v(:,1)= v(:,1).*c(:,1);

LU(:,1)=bc(:,1);     LU(:,2)=h(:,1); LU(:,3)=c(:,1); LU(:,4)=aa(:,1); LU(:,5)=v(:,1);

for j = 1:ny
for k = 1:nz
% Step 2
sum1 = LU(1,2).*RHS(1,j,k);
for i = 2:nx-1
RHS(i,j,k) = RHS(i,j,k) - LU(i,1).*RHS(i-1,j,k);
sum1 = sum1 + LU(i,2)*RHS(i,j,k);
end 
RHS(nx,j,k) = RHS(nx,j,k) - sum1;

%Step 3
RHS(nx,j,k)   = RHS(nx,j,k).*LU(nx,3);

RHS(nx-1,j,k) =  RHS(nx-1,j,k) .* LU(nx-1,3) - RHS(nx,j,k) .* LU(nx-1,5);
for i = nx-2:-1:1
RHS(i,j,k) =  RHS(i,j,k) .* LU(i,3)- RHS(i+1,j,k) .* LU(i,4)- RHS(nx,j,k) .* LU(i,5);
end 
end 
end 


























