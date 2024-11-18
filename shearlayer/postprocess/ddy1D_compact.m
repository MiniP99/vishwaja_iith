function [RHS]=ddy1D_compact(f,dy,nx,ny,nz)

alpha06d1=1/3;a06d1=(14/9)/2;b06d1=(1/9)/4;  %first der coeff
alpha=3;p=17/6;q=3/2;r=3/2;s=-1/6;  %edge nodes
qhat=a06d1;rhat=b06d1;alpha_hat=alpha06d1; %interior nodes
q_p=3/4;alpha_p=1/4;  %node2
alpha_pp = ((40*alpha_hat-1)*q+ 7*(4*alpha_hat-1)*s)/(16*(alpha_hat + 2)*q + 8*(1-4*alpha_hat)*s);
q_pp= (1/3)*(alpha_pp + 2);r_pp= (1/12)*(4*alpha_pp-1);s_pp=0;  %node-3
w1 = (2*alpha_hat + 1)/(2*(q + s));w2 =((8*alpha_hat + 7)*q - 6*(2*alpha_hat + 1)*r+ (8*alpha_hat + 7)*s)/(9*(q + s));
w3 = (4*(alpha_hat + 2)*q + 2*(1 - 4*alpha_hat)*s)/(9*(q + s));     %weights

RHS=zeros(1,ny);

a06 = a06d1/dy;
b06 = b06d1/dy;

a_np_3 = w3*q_pp/dy;
b_np_3 = w3*r_pp/dy;

a_np_2 = w2*q_p/dy;

a_np_1 = w1*( -p/dy);
b_np_1 = w1*(  q/dy);
c_np_1 = w1*(  r/dy);
d_np_1 = w1*(  s/dy);
    RHS(1,1)      =   a_np_1* f(1,1) +  b_np_1*f(1,2) +   c_np_1* f(1,3) +  d_np_1*f(1,4);
    RHS(1,2)      =   a_np_2*(f(1,3) -f(1,1));
    RHS(1,3)      =   a_np_3*(f(1,4) -f(1,2)) + b_np_3*(f(1,5)-f(1,1));
    RHS(1,4:ny-3) =   b06*(f(1,6:ny-1) -f(1,2:ny-5)) + a06*(f(1,5:ny-2) -f(1,3:ny-4));
    RHS(1,ny-2)   =   a_np_3*(f(1,ny-1) -f(1,ny-3)) +b_np_3*(f(1,ny) -f(1,ny-4));
    RHS(1,ny-1)   =   a_np_2*(f(1,ny) -f(1,ny-2));
    RHS(1,ny)     =  -a_np_1* f(1,ny)-b_np_1*f(1,ny-1)-c_np_1* f(1,ny-2)-d_np_1*f(1,ny-3);


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

RHS(1,1) = RHS(1,1).*Tri1(1,2);
for j = 2:ny
    RHS(1,j) = RHS(1,j).*Tri1(j,2) -RHS(1,j-1).*Tri1(j,1);
end
for j = ny-1:-1:1
    RHS(1,j) = RHS(1,j) - Tri1(j,3).*RHS(1,j+1);
end

