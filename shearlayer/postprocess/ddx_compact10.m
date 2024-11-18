function [RHS]=ddx_compact10(f,dx,nx,ny,nz)
alpha10d1=  1.0/2.0;
beta10d1 =  1.0/ 20.0;
a10d1    = (17.0/ 12.0)/ 2.0;
b10d1    = (101.0/150.0)/ 4.0;
c10d1    = (1.0/100.0)/6.0;
a10 = a10d1 * (1/dx);
b10 = b10d1 * (1/dx);
c10 = c10d1 * (1/dx);


RHS = zeros(nx,ny,nz);
     
RHS(1,:,:) = a10 * (f(2,:,:) - f(nx,:,:)) + b10 * (f(3,:,:) - f(nx-1,:,:)) + c10 * (f(4,:,:) - f(nx-2,:,:));    
RHS(2,:,:) = a10 * (f(3,:,:) - f(1,:,:)) + b10 * (f(4,:,:) - f(nx,:,:)) + c10 * (f(5,:,:) - f(nx-1,:,:));
RHS(3,:,:) = a10 * (f(4,:,:) - f(2,:,:)) + b10 * (f(5,:,:) - f(1,:,:)) + c10 * (f(6,:,:) - f(nx,:,:));
RHS(4:nx-3,:,:) = a10 * (f(5:nx-2,:,:) - f(3:nx-4,:,:)) + b10 * (f(6:nx-1,:,:) - f(2:nx-5,:,:)) + c10 * (f(7:nx,:,:) - f(1:nx-6,:,:));
RHS(nx-2,:,:) = a10 * (f(nx-1,:,:) - f(nx-3,:,:)) + b10 * (f(nx,:,:) - f(nx-4,:,:)) + c10 * (f(1,:,:) - f(nx-5,:,:));
RHS(nx-1,:,:) = a10 * (f(nx,:,:) - f(nx-2,:,:)) + b10 * (f(1,:,:) - f(nx-3,:,:)) + c10 * (f(2,:,:) - f(nx-4,:,:));
RHS(nx,:,:) = a10 * (f(1,:,:) - f(nx-1,:,:)) + b10 * (f(2,:,:) - f(nx-2,:,:)) + c10 * (f(3,:,:) - f(nx-3,:,:));

LU = zeros(nx, 9);
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
for i = 3:nx-3
    b(i) = (a - (e / g(i-2)) * h(i-2)) / g(i-1);
    h(i) = c - b(i) * f;
    g(i) = d - (e / g(i-2)) * f - b(i) * h(i-1);
end

% Step 3
b(nx-2) = (a - (e / g(nx-4)) * h(nx-4)) / g(nx-3);
g(nx-2) = d - (e / g(nx-4)) * f - b(nx-2) * h(nx-3);

% Step 4
for i = 3:nx-4
    k(i) = -(k(i-2) * f + k(i-1) * h(i-1)) / g(i);
    v(i) = -(e / g(i-2)) * v(i-2) - b(i) * v(i-1);
end

% Step 5
k(nx-3) = (e - k(nx-5) * f - k(nx-4) * h(nx-4)) / g(nx-3);
k(nx-2) = (a - k(nx-4) * f - k(nx-3) * h(nx-3)) / g(nx-2);
v(nx-3) = f - (e / g(nx-5)) * v(nx-5) - b(nx-3) * v(nx-4);
v(nx-2) = c - (e / g(nx-4)) * v(nx-4) - b(nx-2) * v(nx-3);
g(nx-1) = d - sum(k(1:nx-2) .* v(1:nx-2));

% Step 6
for i = 3:nx-3
    w(i) = -(e / g(i-2)) * w(i-2) - b(i) * w(i-1);
    l(i) = -(l(i-2) * f + l(i-1) * h(i-1)) / g(i);
end

% Step 7
w(nx-2) = f - (e / g(nx-4)) * w(nx-4) - b(nx-2) * w(nx-3);
w(nx-1) = c - sum(k(1:nx-2) .* w(1:nx-2));
l(nx-2) = (e - l(nx-4) * f - l(nx-3) * h(nx-3)) / g(nx-2);
l(nx-1) = (a - sum(l(1:nx-2) .* v(1:nx-2))) / g(nx-1);
g(nx) = d - sum(l(1:nx-1) .* w(1:nx-1));

% Set eg(i) = e/g(i-2)
eg(3:nx-2) = e / g(1:nx-4);

% Set ff = f
ff(1:nx-4) = f;

% Set g = 1/g
g = 1 ./ g;

LU(:, 1) = b;  LU(:, 2) = eg;  LU(:, 3) = k;
LU(:, 4) = l;  LU(:, 5) = g;   LU(:, 6) = h;
LU(:, 7) = ff; LU(:, 8) = v;   LU(:, 9) = w;

for j = 1:ny
    for k = 1:nz
        % Step 8 (update RHS instead of creating z)
        RHS(2, j, k) = RHS(2, j, k) - LU(2, 1) * RHS(1, j, k);
        sum1 = LU(1, 3) * RHS(1, j, k) + LU(2, 3) * RHS(2, j, k);
        sum2 = LU(1, 4) * RHS(1, j, k) + LU(2, 4) * RHS(2, j, k);
        
        % Step 9
        for i = 3:nx-2
            RHS(i, j, k) = RHS(i, j, k) - LU(i, 1) * RHS(i-1, j, k) - LU(i, 2) * RHS(i-2, j, k);
            sum1 = sum1 + LU(i, 3) * RHS(i, j, k);
            sum2 = sum2 + LU(i, 4) * RHS(i, j, k);
        end
        
        % Step 10
        RHS(nx-1, j, k) = RHS(nx-1, j, k) - sum1;
        RHS(nx, j, k) = (RHS(nx, j, k) - sum2 - LU(nx-1, 4) * RHS(nx-1, j, k)) * LU(nx, 5);
        
        % Step 11
        RHS(nx-1, j, k) = (RHS(nx-1, j, k) - LU(nx-1, 9) * RHS(nx, j, k)) * LU(nx-1, 5);
        RHS(nx-2, j, k) = (RHS(nx-2, j, k) - LU(nx-2, 8) * RHS(nx-1, j, k) - LU(nx-2, 9) * RHS(nx, j, k)) * LU(nx-2, 5);
        RHS(nx-3, j, k) = (RHS(nx-3, j, k) - LU(nx-3, 6) * RHS(nx-2, j, k) - LU(nx-3, 8) * RHS(nx-1, j, k) - LU(nx-3, 9) * RHS(nx, j, k)) * LU(nx-3, 5);
        
        for i = nx-4:-1:1
            RHS(i, j, k) = (RHS(i, j, k) - LU(i, 6) * RHS(i+1, j, k) - LU(i, 7) * RHS(i+2, j, k) - LU(i, 8) * RHS(nx-1, j, k) - LU(i, 9) * RHS(nx, j, k)) * LU(i, 5);
        end
    end
end                    
