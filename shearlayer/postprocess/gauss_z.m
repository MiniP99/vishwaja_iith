function [fil] = gauss_z(f,nx,ny,nz)

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

for i = 1:nx
    for j = 1:ny
        fil(i, j, 1) = agf * (f(i, j, 1)) + bgf * (f(i, j, 2) + f(i, j, nz)) ...
            + cgf * (f(i, j, 3) + f(i, j, nz - 1)) + dgf * (f(i, j, 4) + f(i, j, nz - 2)) ...
            + egf * (f(i, j, 5) + f(i, j, nz - 3));
        
        fil(i, j, 2) = agf * (f(i, j, 2)) + bgf * (f(i, j, 3) + f(i, j, 1)) ...
            + cgf * (f(i, j, 4) + f(i, j, nz)) + dgf * (f(i, j, 5) + f(i, j, nz - 1)) ...
            + egf * (f(i, j, 6) + f(i, j, nz - 2));
        
        fil(i, j, 3) = agf * (f(i, j, 3)) + bgf * (f(i, j, 4) + f(i, j, 2)) ...
            + cgf * (f(i, j, 5) + f(i, j, 1)) + dgf * (f(i, j, 6) + f(i, j, nz)) ...
            + egf * (f(i, j, 7) + f(i, j, nz - 1));
        
        fil(i, j, 4) = agf * (f(i, j, 4)) + bgf * (f(i, j, 5) + f(i, j, 3)) ...
            + cgf * (f(i, j, 6) + f(i, j, 2)) + dgf * (f(i, j, 7) + f(i, j, 1)) ...
            + egf * (f(i, j, 8) + f(i, j, nz));
        
        fil(i, j, 5:nz-4) = agf * (f(i,j, 5:nz-4)) + bgf * (f(i,j, 6:nz-3) + f(i, j, 4:nz-5)) ...
            + cgf * (f(i, j, 7:nz-2) + f(i, j, 3:nz-6)) + dgf * (f(i, j, 8:nz-1) + f(i, j, 2:nz-7)) ...
            + egf * (f(i, j, 9:nz) + f(i, j, 1:nz-8));
        
        fil(i, j, nz-3) = agf * (f(i, j, nz-3)) + bgf * (f(i, j, nz-2) + f(i, j, nz-4)) ...
            + cgf * (f(i, j, nz-1) + f(i, j, nz-5)) + dgf * (f(i, j, nz) + f(i, j, nz-6)) ...
            + egf * (f(i, j, 1) + f(i, j, nz-7));
        
        fil(i, j, nz-2) = agf * (f(i, j, nz-2)) + bgf * (f(i, j, nz-1) + f(i, j, nz-3)) ...
            + cgf * (f(i, j, nz) + f(i, j, nz-4)) + dgf * (f(i, j, 1) + f(i, j, nz-5)) ...
            + egf * (f(i, j, 2) + f(i, j, nz-6));
        
        fil(i, j, nz-1) = agf * (f(i, j, nz-1)) + bgf * (f(i, j, nz) + f(i, j, nz-2)) ...
            + cgf * (f(i, j, 1) + f(i,j, nz-3)) + dgf * (f(i, j, 2) + f(i, j, nz-4)) ...
            + egf * (f(i, j, 3) + f(i, j, nz-5));
        
        fil(i, j, nz) = agf * (f(i, j, nz)) + bgf * (f(i, j, 1) + f(i, j, nz-1)) ...
            + cgf * (f(i, j, 2) + f(i, j, nz-2)) + dgf * (f(i, j, 3) + f(i, j, nz-3)) ...
            + egf * (f(i, j, 4) + f(i, j, nz-4));
    end
end
