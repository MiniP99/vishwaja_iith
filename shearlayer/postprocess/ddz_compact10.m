function [RHS]=ddz_compact10(f,dz,nx,ny,nz)
alpha10d1=  1.0/2.0;
beta10d1 =  1.0/ 20.0;
a10d1    = (17.0/ 12.0)/ 2.0;
b10d1    = (101.0/150.0)/ 4.0;
c10d1    = (1.0/100.0)/6.0;
a10 = a10d1 * (1/dz);
b10 = b10d1 * (1/dz);
c10 = c10d1 * (1/dz);

% Assuming n3, n2, and nz are defined

RHS = zeros(nx, ny, nz);
RHS(:,:,1) = a10 * (f(:,:,2) - f(:,:,nz)) + b10 * (f(:,:,3) - f(:,:,nz-1)) + c10 * (f(:,:,4) - f(:,:,nz-2));
RHS(:,:,2) = a10 * (f(:,:,3) - f(:,:,1)) + b10 * (f(:,:,4) - f(:,:,nz)) + c10 * (f(:,:,5) - f(:,:,nz-1));
RHS(:,:,3) = a10 * (f(:,:,4) - f(:,:,2)) + b10 * (f(:,:,5) - f(:,:,1)) + c10 * (f(:,:,6) - f(:,:,nz));
RHS(:,:,4:nz-3) = a10 * (f(:,:,5:nz-2) - f(:,:,3:nz-4)) + b10 * (f(:,:,6:nz-1) - f(:,:,2:nz-5)) + c10 * (f(:,:,7:nz) - f(:,:,1:nz-6));
RHS(:,:,nz-2) = a10 * (f(:,:,nz-1) - f(:,:,nz-3)) + b10 * (f(:,:,nz) - f(:,:,nz-4)) + c10 * (f(:,:,1) - f(:,:,nz-5));
RHS(:,:,nz-1) = a10 * (f(:,:,nz) - f(:,:,nz-2)) + b10 * (f(:,:,1) - f(:,:,nz-3)) + c10 * (f(:,:,2) - f(:,:,nz-4));
RHS(:,:,nz) = a10 * (f(:,:,1) - f(:,:,nz-1)) + b10 * (f(:,:,2) - f(:,:,nz-2)) + c10 * (f(:,:,3) - f(:,:,nz-3));

LU = zeros(nz, 9);
b = LU(:, 1); eg = LU(:, 2); k = LU(:, 3);
l = LU(:, 4); g = LU(:, 5); h = LU(:, 6);
ff = LU(:, 7); v = LU(:, 8); w = LU(:, 9);

d = 1.0;
a = alpha10d1;
c = alpha10d1;
e = beta10d1;
f = beta10d1;

% Step 1
g(1) = d;
b(2) = a / g(1);
h(1) = c;
k(1) = f / g(1);
w(1) = a;
v(1) = e;
l(1) = c / g(1);
g(2) = d - b(2) * h(1);
k(2) = -k(1) * h(1) / g(2);
w(2) = e - b(2) * w(1);
v(2) = -b(2) * v(1);
l(2) = (f - l(1) * h(1)) / g(2);
h(2) = c - b(2) * f;

% Step 2
for i = 3:nz-3
    b(i) = (a - (e / g(i-2)) * h(i-2)) / g(i-1);
    h(i) = c - b(i) * f;
    g(i) = d - (e / g(i-2)) * f - b(i) * h(i-1);
end

% Step 3
b(nz-2) = (a - (e / g(nz-4)) * h(nz-4)) / g(nz-3);
g(nz-2) = d - (e / g(nz-4)) * f - b(nz-2) * h(nz-3);

% Step 4
for i = 3:nz-4
    k(i) = -(k(i-2) * f + k(i-1) * h(i-1)) / g(i);
    v(i) = -(e / g(i-2)) * v(i-2) - b(i) * v(i-1);
end

% Step 5
k(nz-3) = (e - k(nz-5) * f - k(nz-4) * h(nz-4)) / g(nz-3);
k(nz-2) = (a - k(nz-4) * f - k(nz-3) * h(nz-3)) / g(nz-2);
v(nz-3) = f - (e / g(nz-5)) * v(nz-5) - b(nz-3) * v(nz-4);
v(nz-2) = c - (e / g(nz-4)) * v(nz-4) - b(nz-2) * v(nz-3);
g(nz-1) = d - sum(k(1:nz-2) .* v(1:nz-2));

% Step 6
for i = 3:nz-3
    w(i) = -(e / g(i-2)) * w(i-2) - b(i) * w(i-1);
    l(i) = -(l(i-2) * f + l(i-1) * h(i-1)) / g(i);
end

% Step 7
w(nz-2) = f - (e / g(nz-4)) * w(nz-4) - b(nz-2) * w(nz-3);
w(nz-1) = c - sum(k(1:nz-2) .* w(1:nz-2));
l(nz-2) = (e - l(nz-4) * f - l(nz-3) * h(nz-3)) / g(nz-2);
l(nz-1) = (a - sum(l(1:nz-2) .* v(1:nz-2))) / g(nz-1);
g(nz) = d - sum(l(1:nz-1) .* w(1:nz-1));

% Set eg(k) = e/g(k-2)
eg(3:nz-2) = e / g(1:nz-4);

% Set ff = f
ff(1:nz-4) = f;

% Set g = 1/g
g = 1 ./ g;

LU(:, 1) = b;  LU(:, 2) = eg;  LU(:, 3) = k;
LU(:, 4) = l;  LU(:, 5) = g;   LU(:, 6) = h;
LU(:, 7) = ff; LU(:, 8) = v;   LU(:, 9) = w;

for i = 1:nx
    for j = 1:ny
        % Step 8 (update RHS instead of creating z)
        RHS(i, j, 2) = RHS(i, j, 2) - LU(2, 1) * RHS(i, j, 1);
        sum1 = LU(1, 3) * RHS(i, j, 1) + LU(2, 3) * RHS(i, j, 2);
        sum2 = LU(1, 4) * RHS(i, j, 1) + LU(2, 4) * RHS(i, j, 2);

        % Step 9
        for k = 3:nz-2
            RHS(i, j, k) = RHS(i, j, k) - LU(k, 1) * RHS(i, j, k-1) - LU(k, 2) * RHS(i, j, k-2);
            sum1 = sum1 + LU(k, 3) * RHS(i, j, k);
            sum2 = sum2 + LU(k, 4) * RHS(i, j, k);
        end

        % Step 10
        RHS(i, j, nz-1) = RHS(i, j, nz-1) - sum1;
        RHS(i, j, nz) = (RHS(i, j, nz) - sum2 - LU(nz-1, 4) * RHS(i, j, nz-1)) * LU(nz, 5);

        % Step 11
        RHS(i, j, nz-1) = (RHS(i, j, nz-1) - LU(nz-1, 9) * RHS(i, j, nz)) * LU(nz-1, 5);
        RHS(i, j, nz-2) = (RHS(i, j, nz-2) - LU(nz-2, 8) * RHS(i, j, nz-1) - LU(nz-2, 9) * RHS(i, j, nz)) * LU(nz-2, 5);
        RHS(i, j, nz-3) = (RHS(i, j, nz-3) - LU(nz-3, 6) * RHS(i, j, nz-2) - LU(nz-3, 8) * RHS(i, j, nz-1) - LU(nz-3, 9) * RHS(i, j, nz)) * LU(nz-3, 5);

        for k = nz-4:-1:1
            RHS(i, j, k) = (RHS(i, j, k) - LU(k, 6) * RHS(i, j, k+1) - LU(k, 7) * RHS(i, j, k+2) - LU(k, 8) * RHS(i, j, nz-1) - LU(k, 9) * RHS(i, j, nz)) * LU(k, 5);
        end
    end
end
