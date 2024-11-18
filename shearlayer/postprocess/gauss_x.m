function [fil] = gauss_x(f,nx,ny,nz)

agf = 3565/10368 ; bgf = 3091/12960 ; cgf = 1997/25920; dgf = 149/12960; egf = 107/103680;
% 1st point
b1_agf = 5/6; b1_bgf = 1/6;
% 2nd point
b2_agf = 2/3 ; b2_bgf = 1/6;
% 3rd point
b3_agf = 31/64 ; b3_bgf = 7/32; b3_cgf = 5/128;
% 4rd point
b4_agf = 17/48 ; b4_bgf = 15/64; b4_cgf = 7/96; b4_dgf = 1/64;

% Periodic boundary for now

% Initialize fil array
fil = zeros(nx, ny, nz);

for j = 1:ny
    for k = 1:nz
        fil(1, j, k) = agf * (f(1, j, k)) + bgf * (f(2, j, k) + f(nx, j, k)) ...
            + cgf * (f(3, j, k) + f(nx - 1, j, k)) + dgf * (f(4, j, k) + f(nx - 2, j, k)) ...
            + egf * (f(5, j, k) + f(nx - 3, j, k));
        
        fil(2, j, k) = agf * (f(2, j, k)) + bgf * (f(3, j, k) + f(1, j, k)) ...
            + cgf * (f(4, j, k) + f(nx, j, k)) + dgf * (f(5, j, k) + f(nx - 1, j, k)) ...
            + egf * (f(6, j, k) + f(nx - 2, j, k));
        
        fil(3, j, k) = agf * (f(3, j, k)) + bgf * (f(4, j, k) + f(2, j, k)) ...
            + cgf * (f(5, j, k) + f(1, j, k)) + dgf * (f(6, j, k) + f(nx, j, k)) ...
            + egf * (f(7, j, k) + f(nx - 1, j, k));
        
        fil(4, j, k) = agf * (f(4, j, k)) + bgf * (f(5, j, k) + f(3, j, k)) ...
            + cgf * (f(6, j, k) + f(2, j, k)) + dgf * (f(7, j, k) + f(1, j, k)) ...
            + egf * (f(8, j, k) + f(nx, j, k));
        
        fil(5:nx-4, j, k) = agf * (f(5:nx-4, j, k)) + bgf * (f(6:nx-3, j, k) + f(4:nx-5, j, k)) ...
            + cgf * (f(7:nx-2, j, k) + f(3:nx-6, j, k)) + dgf * (f(8:nx-1, j, k) + f(2:nx-7, j, k)) ...
            + egf * (f(9:nx, j, k) + f(1:nx-8, j, k));
        
        fil(nx-3, j, k) = agf * (f(nx-3, j, k)) + bgf * (f(nx-2, j, k) + f(nx-4, j, k)) ...
            + cgf * (f(nx-1, j, k) + f(nx-5, j, k)) + dgf * (f(nx, j, k) + f(nx-6, j, k)) ...
            + egf * (f(1, j, k) + f(nx-7, j, k));
        
        fil(nx-2, j, k) = agf * (f(nx-2, j, k)) + bgf * (f(nx-1, j, k) + f(nx-3, j, k)) ...
            + cgf * (f(nx, j, k) + f(nx-4, j, k)) + dgf * (f(1, j, k) + f(nx-5, j, k)) ...
            + egf * (f(2, j, k) + f(nx-6, j, k));
        
        fil(nx-1, j, k) = agf * (f(nx-1, j, k)) + bgf * (f(nx, j, k) + f(nx-2, j, k)) ...
            + cgf * (f(1, j, k) + f(nx-3, j, k)) + dgf * (f(2, j, k) + f(nx-4, j, k)) ...
            + egf * (f(3, j, k) + f(nx-5, j, k));
        
        fil(nx, j, k) = agf * (f(nx, j, k)) + bgf * (f(1, j, k) + f(nx-1, j, k)) ...
            + cgf * (f(2, j, k) + f(nx-2, j, k)) + dgf * (f(3, j, k) + f(nx-3, j, k)) ...
            + egf * (f(4, j, k) + f(nx-4, j, k));
    end
end
