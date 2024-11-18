function [RHS]=ddeta_compact10_gen(f,deta)  %General boundary condition

[nx,ny,nz] = size(f);
alpha10d1=  1.0/2.0; beta10d1 =  1.0/ 20.0;
a10d1    = (17.0/ 12.0)/ 2.0; b10d1    = (101.0/150.0)/ 4.0; c10d1    = (1.0/100.0)/6.0;
a10 = a10d1 * (1/deta); b10 = b10d1 * (1/deta); c10 = c10d1 * (1/deta);

alpha = 3 ; p = -17/6 ; q = 3/2 ; r = 3/2 ; s = -1/6 ; 

q_hat = a10d1; r_hat = b10d1; s_hat = c10d1;
alpha_hat = alpha10d1; beta_hat = beta10d1;

q_p             = 3/4;  alpha_p         = 1/4;

alpha_ppp       = (8*r_hat - 175*s_hat)/(18*r_hat - 550*s_hat); beta_ppp        = (1/20)*(-3 + 8*alpha_ppp);
q_ppp           = (1/12)*(12 - 7*alpha_ppp) ; r_ppp = (1/600)*(568*alpha_ppp - 183) ; s_ppp    = (1/300)*(9*alpha_ppp - 4) ;

alpha_pp        = ((17*(s*(r_hat + 2*s_hat) - q*(q_hat + r_hat+ s_hat)))/(72*(q + s)*(q_hat + r_hat - s_hat*(q_ppp/s_ppp- 1))) - 8/9)/((19*(s*(r_hat + 2*s_hat) - q*(q_hat+ r_hat + s_hat)))/(24*(q + s)*(q_hat + r_hat - s_hat*(q_ppp/s_ppp- 1))) - 1/3);

beta_pp         = (1/12)*(-1 + 3*alpha_pp);   q_pp = (2/18)*(8 - 3*alpha_pp) ;
r_pp            = (1/72)*(-17 + 57*alpha_pp); s_pp  = 0;

w1              = (q_hat + 2*r_hat + 3*s_hat)/(q + s);
w2              = (1/q_p)*(r_hat + s_hat*(1 + q_ppp/s_ppp) - r*(q_hat+ 2*r_hat + 3*s_hat)/(q + s) );
w3              = (q_hat + r_hat + s_hat*(1 - q_ppp/s_ppp))/(r_pp) ;
w4              = s_hat/s_ppp ;

a_np_4 = w4*q_ppp * (1/deta) ; b_np_4 = w4*r_ppp * (1/deta); c_np_4 = w4*s_ppp * (1/deta);
a_np_3 = w3*q_pp * (1/deta);   b_np_3 = w3*r_pp * (1/deta);
a_np_2 = w2*q_p * (1/deta);            
a_np_1 = w1*( p * (1/deta));  b_np_1 = w1*( q * (1/deta)); c_np_1 = w1*( r * (1/deta)); d_np_1 = w1*( s * (1/deta));


RHS = zeros(nx, ny, nz);
RHS(:,1,:) =   a_np_1* f(:,1,:) +  b_np_1*f(:,2,:) +   c_np_1* f(:,3,:) +  d_np_1*f(:,4,:);
RHS(:,2,:) =   a_np_2*(f(:,3,:) - f(:,1,:));
RHS(:,3,:) =   a_np_3*(f(:,4,:) - f(:,2,:)) + b_np_3*(f(:,5,:) - f(:,1,:));
RHS(:,4,:) =   a_np_4*(f(:,5,:) - f(:,3,:)) + b_np_4*(f(:,6,:) - f(:,2,:))+ c_np_4*(f(:,7,:)-f(:,1,:));
RHS(:,5:ny-4,:) = a10 * (f(:,6:ny-3,:) - f(:,4:ny-5,:)) + b10 * (f(:,7:ny-2,:) - f(:,3:ny-6,:)) + c10 * (f(:,8:ny-1,:) - f(:,2:ny-7,:));
RHS(:,ny-3, :) = a_np_4 * (f(:,ny-2,:) - f(:,ny-4,:)) + b_np_4 * (f(:,ny-1,:) - f(:,ny-5,:)) + c_np_4 * (f(:,ny,:) - f(:,ny-6,:));
RHS(:,ny-2, :) = a_np_3 * (f(:,ny-1,:) - f(:,ny-3,:)) + b_np_3 * (f(:,ny,:) - f(:,ny-4,:));
RHS(:,ny-1, :) = a_np_2 * (f(:,ny,:) - f(:,ny-2,:));
RHS(:,ny,:)   = -a_np_1* f(:,ny,:) -  b_np_1*f(:,ny-1,:) - c_np_1* f(:,ny-2,:) -  d_np_1*f(:,ny-3,:);

penta1 = zeros(ny,11);
bt   = penta1(:, 1);    b    = penta1(:, 2);    d    = penta1(:, 3);    a    = penta1(:, 4);
at   = penta1(:, 5);    e    = penta1(:, 6);    obc  = penta1(:, 7);    f    = penta1(:, 8);
g    = penta1(:, 9);    eobc = penta1(:, 10);

at(:) = beta_hat; bt(:) = beta_hat; a(:)  = alpha_hat; b(:)  = alpha_hat; d(:)  = 1.0;

bt(1) = w1*0; b (1) = w1*0; d (1) = w1*1.0; a (1) = w1*alpha; at(1) = w1*0;
bt(2) = w2*0; b (2) = w2*alpha_p; d (2) = w2*1; a (2) = w2*alpha_p; at(2) = w2*0;
bt(3) = w3*beta_pp; b (3) = w3*alpha_pp; d (3) = w3*1; a (3) = w3*alpha_pp; at(3) = w3*beta_pp;
bt(4) = w4*beta_ppp; b (4) = w4*alpha_ppp; d (4) = w4*1; a (4) = w4*alpha_ppp; at(4) = w4*beta_ppp;

bt(ny) = w1*0; b(ny) = w1*alpha; d(ny) = w1*1; a(ny) = w1*0; at(ny) = w1*0;
bt(ny-1) = w2*0; b(ny-1) = w2*alpha_p; d(ny-1) = w2*1; a(ny-1) = w2*alpha_p; at(ny-1) = w2*0;
bt(ny-2) = w3*beta_pp; b(ny-2) = w3*alpha_pp; d(ny-2) = w3*1; a(ny-2) = w3*alpha_pp; at(ny-2) = w3*beta_pp;
bt(ny-3) = w4*beta_ppp; b(ny-3) = w4*alpha_ppp; d(ny-3) = w4*1; a(ny-3) = w4*alpha_ppp; at(ny-3) = w4*beta_ppp;

% Step 1
obc(1) = 1/d(1);

% Step 2
obc(2) = 1/(d(2) - b(2)*a(1)*obc(1));

% Step 3
e(1) = a(1);
f(2) = b(2)*obc(1);
for j = 3:ny
    g(j) = bt(j)*obc(j-2);
    e(j-1) = a(j-1) - f(j-1)*at(j-2);
    f(j) = (b(j) - g(j)*e(j-2))*obc(j-1);
    obc(j) = 1/(d(j) - f(j)*e(j-1) - g(j)*at(j-2));
end

eobc = e.*obc;

penta1(:, 1) = bt;  penta1(:, 2) = b;   penta1(:, 3) = d;   penta1(:, 4) = a;
penta1(:, 5) = at;  penta1(:, 6) = e;   penta1(:, 7) = obc; penta1(:, 8) = f;
penta1(:, 9) = g;   penta1(:, 10) = eobc;

for i = 1:nx
    for k = 1:nz
        % Step 1
        RHS(i, 2, k) = RHS(i, 2, k) - penta1(2, 8) * RHS(i, 1, k);
        for j = 3:ny
            RHS(i, j, k) = RHS(i, j, k) - penta1(j, 9) * RHS(i, j-2, k) - penta1(j, 8) * RHS(i, j-1, k);
        end

        % Step 2
        RHS(i, ny, k) = RHS(i, ny, k) * penta1(ny, 7);

        RHS(i, ny-1, k) = RHS(i, ny-1, k) * penta1(ny-1, 7) - penta1(ny-1, 10) * RHS(i, ny, k);
        for j = ny-2:-1:1
            RHS(i, j, k) = RHS(i, j, k) * penta1(j, 7) - RHS(i, j+2, k) * penta1(j, 5) * penta1(j, 7) - RHS(i, j+1, k) * penta1(j, 10);
        end
    end
end














