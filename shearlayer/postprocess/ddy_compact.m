function [RHS]=ddy_compact(f,dy,nx,ny,nz)
% General boundary condition

alpha06d1=1/3;a06d1=(14/9)/2;b06d1=(1/9)/4;  %first der coeff
alpha=3;p=17/6;q=3/2;r=3/2;s=-1/6;  %edge nodes
qhat=a06d1;rhat=b06d1;alpha_hat=alpha06d1; %interior nodes
q_p=3/4;alpha_p=1/4;  %node2
alpha_pp = ((40*alpha_hat-1)*q+ 7*(4*alpha_hat-1)*s)/(16*(alpha_hat + 2)*q + 8*(1-4*alpha_hat)*s);
q_pp= (1/3)*(alpha_pp + 2);r_pp= (1/12)*(4*alpha_pp-1);s_pp=0;  %node-3
w1 = (2*alpha_hat + 1)/(2*(q + s));w2 =((8*alpha_hat + 7)*q - 6*(2*alpha_hat + 1)*r+ (8*alpha_hat + 7)*s)/(9*(q + s));
w3 = (4*(alpha_hat + 2)*q + 2*(1 - 4*alpha_hat)*s)/(9*(q + s));     %weights

RHS=zeros(nx,ny,nz);

a06 = a06d1/dy;
b06 = b06d1/dy;

a_np_3 = w3*q_pp/dy;
b_np_3 = w3*r_pp/dy;

a_np_2 = w2*q_p/dy;

a_np_1 = w1*( -p/dy);
b_np_1 = w1*(  q/dy);
c_np_1 = w1*(  r/dy);
d_np_1 = w1*(  s/dy);
for k = 1:nz
    RHS(:,1,k)      =   a_np_1* f(:,1,k) +  b_np_1*f(:,2,k) +   c_np_1* f(:,3,k) +  d_np_1*f(:,4,k);
    RHS(:,2,k)      =   a_np_2*(f(:,3,k) -f(:,1,k));
    RHS(:,3,k)      =   a_np_3*(f(:,4,k) -f(:,2,k)) + b_np_3*(f(:,5,k)-f(:,1,k));
    RHS(:,4:ny-3,k) =   b06*(f(:,6:ny-1,k) -f(:,2:ny-5,k)) + a06*(f(:,5:ny-2,k) -f(:,3:ny-4,k));
    RHS(:,ny-2  ,k) =   a_np_3*(f(:,ny-1,k) -f(:,ny-3,k)) +b_np_3*(f(:,ny,k) -f(:,ny-4,k));
    RHS(:,ny-1  ,k) =   a_np_2*(f(:,ny,k) -f(:,ny-2,k));
    RHS(:,ny    ,k) =  -a_np_1* f(:,ny,k)-b_np_1*f(:,ny-1,k)-c_np_1* f(:,ny-2,k)-d_np_1*f(:,ny-3,k);
end

a  = ones(1,ny)*alpha_hat;   b=ones(1,ny);      c  = ones(1,ny)*alpha_hat;
a (1,1) = w1*0;
a (1,2) = w2*alpha_p;
a (1,3) = w3*alpha_pp;

b (1,1) = w1;
b (1,2) = w2*1;
b (1,3) = w3*1;

c (1,1) = w1*alpha;
c (1,2) = w2*alpha_p;
c (1,3) = w3*alpha_pp;

c (1,ny  ) = w1*0;
c (1,ny-1) = w2*alpha_p;
c (1,ny-2) = w3*alpha_pp;

b (1,ny  ) = w1*1;
b (1,ny-1) = w2*1;
b (1,ny-2) = w3*1;

a (1,ny  ) = w1*alpha;
a (1,ny-1) = w2*alpha_p;
a (1,ny-2) = w3*alpha_pp;

cp =zeros(1,ny);
den=zeros(1,ny);
cp(1,1) = c(1,1)/b(1,1);
for j = 2:ny-1
    cp(1,j) = c(1,j)/(b(1,j) - a(1,j)*cp(1,j-1));
end

den(1,1) = 1/b(1,1);
den(1,2:ny) = 1./(b(1,2:ny) - a(1,2:ny).*cp(1,1:ny-1));

Tri1=zeros(ny,3);
Tri1(:,1) = a(1,:).*den(1,:);
Tri1(:,2) = den(1,:);
Tri1(:,3) = cp(1,:);

for i = 1:nx
    for k = 1:nz
        RHS(i,1,k) = RHS(i,1,k).*Tri1(1,2);
        for j = 2:ny
            RHS(i,j,k) = RHS(i,j,k).*Tri1(j,2) -RHS(i,j-1,k).*Tri1(j,1);
        end
        for j = ny-1:-1:1
            RHS(i,j,k) = RHS(i,j,k) - Tri1(j,3).*RHS(i,j+1,k);
        end
    end
end


