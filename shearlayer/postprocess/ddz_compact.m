function [RHS]=ddz_compact(f,dz,nx,ny,nz)

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

%Periodic for now



a06 = a06d1/dz;
b06 = b06d1/dz;

RHS=zeros(nx,ny,nz);
RHS(:,:,1) =      a06 *(f(:,:,2)- f(:,:,nz))        + b06*(f(:,:,3)-f(:,:,nz-1));
RHS(:,:,2) =      a06 *(f(:,:,3)- f(:,:,1))         + b06*(f(:,:,4)-f(:,:,nz));
RHS(:,:,3:nz-2) = a06*(f(:,:,4:nz-1)-f(:,:,2:nz-3)) + b06*(f(:,:,5:nz)-f(:,:,1:nz-4));
RHS(:,:,nz-1) =   a06*( f(:,:,nz)-f(:,:,nz-2))      + b06*(f(:,:,1)-f(:,:,nz-3));
RHS(:,:,nz)   =   a06 *(f(:,:,1)- f(:,:,nz-1))      + b06*( f(:,:,2)- f(:,:,nz-2) );


LU = zeros(nz,5);
bc=LU(:,1);     h=LU(:,2); c=LU(:,3); aa=LU(:,4); v=LU(:,5);
b=alpha06d1;d=1;a=alpha06d1;

c(1,1) = d;
v(1,1) = b;
h(1,1) = a/c(1,1);

% Step 1
for k = 2:nz-1
c(k,1) = d - ( b/c(k-1,1) ).*a;
end 
for k = 2:nz-2
v(k,1) = -( b/c(k-1,1) ).*v(k-1,1);
h(k,1) = -( a/c(k,1) ).*h(k-1,1);
end 
v(nz-1,1) = a - ( b/c(nz-2,1) ).*v(nz-2,1);
h(nz-1,1) = ( b - h(nz-2,1).*a ) ./ c(nz-1,1);
c(nz,1)   = d - sum( h(1:nz-1,1).*v(1:nz-1,1) );

bc(2:nz-1,1) = b ./ c(1:nz-2,1);
aa(1:nz-2,1) = a;
c(:,1) = 1./c(:,1);

% Overwrite aa by aa*c
aa(:,1) = aa(:,1).*c(:,1);

%Overwrite v by v*c
v(:,1)= v(:,1).*c(:,1);

LU(:,1)=bc(:,1);     LU(:,2)=h(:,1); LU(:,3)=c(:,1); LU(:,4)=aa(:,1); LU(:,5)=v(:,1);

for i = 1:nx
for j = 1:ny
% Step 2
sum1 = LU(1,2).*RHS(i,j,1);
for k = 2:nz-1
RHS(i,j,k) = RHS(i,j,k) - LU(k,1).*RHS(i,j,k-1);
sum1 = sum1 + LU(k,2)*RHS(i,j,k);
end 
RHS(i,j,nz) = RHS(i,j,nz) - sum1;
%Step 3
RHS(i,j,nz)   = RHS(i,j,nz).*LU(nz,3);
RHS(i,j,nz-1) =  RHS(i,j,nz-1) .* LU(nz-1,3) - RHS(i,j,nz) .* LU(nz-1,5);
for k = nz-2:-1:1
RHS(i,j,k) =  RHS(i,j,k) .* LU(k,3)- RHS(i,j,k+1) .* LU(k,4)- RHS(i,j,nz) .* LU(k,5);
end 
end 
end 
