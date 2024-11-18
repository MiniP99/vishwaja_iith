function [RHS]=ddzeta_compact10_sym(f,dzeta)  %symmetric boundary condition
[nx,ny,nz]=size(f);
alpha10d1=  1.0/2.0;
beta10d1 =  1.0/ 20.0;
a10d1    = (17.0/ 12.0)/ 2.0;
b10d1    = (101.0/150.0)/ 4.0;
c10d1    = (1.0/100.0)/6.0;
a10 = a10d1 * (1/dzeta);
b10 = b10d1 * (1/dzeta);
c10 = c10d1 * (1/dzeta);

q_hat = a10d1;
r_hat = b10d1;
s_hat = c10d1;
alpha_hat = alpha10d1;
beta_hat = beta10d1;


RHS = zeros(nx, ny, nz);
RHS(:, :, 1) = 0;
RHS(:, :, 2) = a10 * (f(:, :, 3) - f(:, :, 1)) + b10 * (f(:, :, 4) - f(:, :, 2)) + c10 * (f(:, :, 5) - f(:, :, 3));
RHS(:, :, 3) = a10 * (f(:, :, 4) - f(:, :, 2)) + b10 * (f(:, :, 5) - f(:, :, 1)) + c10 * (f(:, :, 6) - f(:, :, 2));
RHS(:, :, 4) = a10 * (f(:, :, 5) - f(:, :, 3)) + b10 * (f(:, :, 6) - f(:, :, 2)) + c10 * (f(:, :, 7) - f(:, :, 1));
RHS(:, :, 5:nz-4) = a10 * (f(:, :, 6:nz-3) - f(:, :, 4:nz-5)) + b10 * (f(:, :, 7:nz-2) - f(:, :, 3:nz-6)) + c10 * (f(:, :, 8:nz-1) - f(:, :, 2:nz-7));
RHS(:, :, nz-3) = a10 * (f(:, :, nz-2) - f(:, :, nz-4)) + b10 * (f(:, :, nz-1) - f(:, :, nz-5)) + c10 * (f(:, :, nz) - f(:, :, nz-6));
RHS(:, :, nz-2) = a10 * (f(:, :, nz-1) - f(:, :, nz-3)) + b10 * (f(:, :, nz) - f(:, :, nz-4)) + c10 * (f(:, :, nz-1) - f(:, :, nz-5));
RHS(:, :, nz-1) = a10 * (f(:, :, nz) - f(:, :, nz-2)) + b10 * (f(:, :, nz-1) - f(:, :, nz-3)) + c10 * (f(:, :, nz-2) - f(:, :, nz-4));
RHS(:, :, nz) = 0;


penta1 = zeros(nz,11);
bt   = penta1(:, 1);    b    = penta1(:, 2);    d    = penta1(:, 3);    a    = penta1(:, 4);
at   = penta1(:, 5);    e    = penta1(:, 6);    obc  = penta1(:, 7);    f    = penta1(:, 8);
g    = penta1(:, 9);    eobc = penta1(:, 10);

at(:) = beta_hat;
bt(:) = beta_hat;
a(:)  = alpha_hat;
b(:)  = alpha_hat;
d(:)  = 1.0;

bt(1) = 0;      b(1) = 0;         d(1) = 1;            a(1) = 0;         at(1) = 0;
bt(2) = 0;      b(2) = alpha_hat; d(2) = 1 - beta_hat; a(2) = alpha_hat; at(2) = beta_hat;
bt(nz) = 0;     b(nz) = 0;        d(nz) = 1;           a(nz) = 0;        at(nz) = 0;
bt(nz-1) = beta_hat; b(nz-1) = alpha_hat; d(nz-1) = 1 - beta_hat; a(nz-1) = alpha_hat; at(nz-1) = 0;

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
