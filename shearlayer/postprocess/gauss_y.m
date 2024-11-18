function [fil] = gauss_y(f,nx,ny,nz)

agf = 3565/10368 ; bgf = 3091/12960 ; cgf = 1997/25920; dgf = 149/12960; egf = 107/103680;
% 1st point
b1_agf = 5/6; b1_bgf = 1/6;
% 2nd point
b2_agf = 2/3 ; b2_bgf = 1/6;
% 3rd point
b3_agf = 31/64 ; b3_bgf = 7/32; b3_cgf = 5/128;
% 4rd point
b4_agf = 17/48 ; b4_bgf = 15/64; b4_cgf = 7/96; b4_dgf = 1/64;

% Initialize fil array
fil = zeros(nx, ny, nz);

% General boundary condition
for k = 1:nz
    fil(:, 1, k) = b1_agf * (f(:, 1, k)) + b1_bgf * (f(:, 2, k));
    fil(:, 2, k) = b2_agf * (f(:, 2, k)) + b2_bgf * (f(:, 3, k) + f(:, 1, k));
    fil(:, 3, k) = b3_agf * (f(:, 3, k)) + b3_bgf * (f(:, 4, k) + f(:, 2, k)) + b3_cgf * (f(:, 5, k) + f(:, 1, k));
    fil(:, 4, k) = b4_agf * (f(:, 4, k)) + b4_bgf * (f(:, 5, k) + f(:, 3, k)) + ...
                   b4_cgf * (f(:, 6, k) + f(:, 2, k)) + b4_dgf * (f(:, 7, k) + f(:, 1, k));
    fil(:, 5:ny-4, k) = agf * (f(:, 5:ny-4, k)) + bgf * (f(:, 6:ny-3, k) + f(:, 4:ny-5, k)) ...
        + cgf * (f(:, 7:ny-2, k) + f(:, 3:ny-6, k)) + dgf * (f(:, 8:ny-1, k) + f(:, 2:ny-7, k)) ...
        + egf * (f(:, 9:ny, k) + f(:, 1:ny-8, k));

    fil(:, ny-3, k) = b4_agf * (f(:, ny-3, k)) + b4_bgf * (f(:, ny-2, k) + f(:, ny-4, k)) ...
        + b4_cgf * (f(:, ny-1, k) + f(:, ny-5, k)) + b4_dgf * (f(:, ny, k) + f(:, ny-6, k));

    fil(:, ny-2, k) = b3_agf * (f(:, ny-2, k)) + b3_bgf * (f(:, ny-1, k) + f(:, ny-3, k)) ...
        + b3_cgf * (f(:, ny, k) + f(:, ny-4, k));

    fil(:, ny-1, k) = b2_agf * (f(:, ny-1, k)) + b2_bgf * (f(:, ny, k) + f(:, ny-2, k));

    fil(:, ny, k) = b1_agf * (f(:, ny, k)) + b1_bgf * (f(:, ny-1, k));
end
