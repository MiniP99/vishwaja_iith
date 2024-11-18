function [RHS]=ddzeta_compact10_gen(f,dzeta)  %General boundary condition
[nx,ny,nz]=size(f);
alpha10d1=  1.0/2.0; beta10d1 =  1.0/ 20.0;
a10d1    = (17.0/ 12.0)/ 2.0; b10d1    = (101.0/150.0)/ 4.0; c10d1    = (1.0/100.0)/6.0;
a10 = a10d1 * (1/dzeta); b10 = b10d1 * (1/dzeta); c10 = c10d1 * (1/dzeta);

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

a_np_4 = w4*q_ppp * (1/dzeta) ; b_np_4 = w4*r_ppp * (1/dzeta); c_np_4 = w4*s_ppp * (1/dzeta);
a_np_3 = w3*q_pp * (1/dzeta);   b_np_3 = w3*r_pp * (1/dzeta);
a_np_2 = w2*q_p * (1/dzeta);            
a_np_1 = w1*( p * (1/dzeta));  b_np_1 = w1*( q * (1/dzeta)); c_np_1 = w1*( r * (1/dzeta)); d_np_1 = w1*( s * (1/dzeta));


RHS = zeros(nx, ny, nz);
RHS(:,:,1) =   a_np_1* f(:,:,1) +  b_np_1*f(:,:,2) +   c_np_1* f(:,:,3) +  d_np_1*f(:,:,4);
RHS(:,:,2) =   a_np_2*(f(:,:,3) - f(:,:,1));
RHS(:,:,3) =   a_np_3*(f(:,:,4) - f(:,:,2)) + b_np_3*(f(:,:,5) - f(:,:,1));
RHS(:,:,4) =   a_np_4*(f(:,:,5) - f(:,:,3)) + b_np_4*(f(:,:,6) - f(:,:,2))+ c_np_4*(f(:,:,7)-f(:,:,1));
RHS(:,:,5:nz-4) = a10 * (f(:,:,6:nz-3) - f(:,:,4:nz-5)) + b10 * (f(:,:,7:nz-2) - f(:,:,3:nz-6)) + c10 * (f(:,:,8:nz-1) - f(:,:,2:nz-7));
RHS(:,:,nz-3) = a_np_4 * (f(:,:,nz-2) - f(:,:,nz-4)) + b_np_4 * (f(:,:,nz-1) - f(:,:,nz-5)) + c_np_4 * (f(:,:,nz) - f(:,:,nz-6));
RHS(:,:,nz-2) = a_np_3 * (f(:,:,nz-1) - f(:,:,nz-3)) + b_np_3 * (f(:,:,nz) - f(:,:,nz-4));
RHS(:,:,nz-1) = a_np_2 * (f(:,:,nz) - f(:,:,nz-2));
RHS(:,:,nz)   = -a_np_1* f(:,:,nz) -  b_np_1*f(:,:,nz-1) - c_np_1* f(:,:,nz-2) -  d_np_1*f(:,:,nz-3);

penta1 = zeros(ny,11);
bt   = penta1(:, 1);    b    = penta1(:, 2);    d    = penta1(:, 3);    a    = penta1(:, 4);
at   = penta1(:, 5);    e    = penta1(:, 6);    obc  = penta1(:, 7);    f    = penta1(:, 8);
g    = penta1(:, 9);    eobc = penta1(:, 10);

at(:) = beta_hat; bt(:) = beta_hat; a(:)  = alpha_hat; b(:)  = alpha_hat; d(:)  = 1.0;

bt(1) = w1*0; b (1) = w1*0; d (1) = w1*1.0; a (1) = w1*alpha; at(1) = w1*0;
bt(2) = w2*0; b (2) = w2*alpha_p; d (2) = w2*1; a (2) = w2*alpha_p; at(2) = w2*0;
bt(3) = w3*beta_pp; b (3) = w3*alpha_pp; d (3) = w3*1; a (3) = w3*alpha_pp; at(3) = w3*beta_pp;
bt(4) = w4*beta_ppp; b (4) = w4*alpha_ppp; d (4) = w4*1; a (4) = w4*alpha_ppp; at(4) = w4*beta_ppp;

bt(nz) = w1*0; b(nz) = w1*alpha; d(nz) = w1*1; a(nz) = w1*0; at(nz) = w1*0;
bt(nz-1) = w2*0; b(nz-1) = w2*alpha_p; d(nz-1) = w2*1; a(nz-1) = w2*alpha_p; at(nz-1) = w2*0;
bt(nz-2) = w3*beta_pp; b(nz-2) = w3*alpha_pp; d(nz-2) = w3*1; a(nz-2) = w3*alpha_pp; at(nz-2) = w3*beta_pp;
bt(nz-3) = w4*beta_ppp; b(nz-3) = w4*alpha_ppp; d(nz-3) = w4*1; a(nz-3) = w4*alpha_ppp; at(nz-3) = w4*beta_ppp;

% Step 1
obc(1) = 1/d(1);

% Step 2
obc(2) = 1/(d(2) - b(2)*a(1)*obc(1));

% Step 3
e(1) = a(1);
f(2) = b(2)*obc(1);
for j = 3:nz
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
    for j = 1:ny
        % Step 1
        RHS(i, j, 2) = RHS(i, j, 2) - penta1(2, 8) * RHS(i, j, 1);
        for k = 3:nz
            RHS(i, j, k) = RHS(i, j, k) - penta1(k, 9) * RHS(i, j, k-2) - penta1(k, 8) * RHS(i, j, k-1);
        end

        % Step 2
        RHS(i, j, nz) = RHS(i, j, nz) * penta1(nz, 7);

        RHS(i, j, nz-1) = RHS(i, j, nz-1) * penta1(nz-1, 7) - penta1(nz-1, 10) * RHS(i, j, nz);
        for k = nz-2:-1:1
            RHS(i, j, k) = RHS(i, j, k) * penta1(k, 7) - RHS(i, j, k+2) * penta1(k, 5) * penta1(k, 7) - RHS(i, j, k+1) * penta1(k, 10);
        end
    end
end

