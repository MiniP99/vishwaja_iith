function [k1,k2,k3,kxdeal,kydeal,kzdeal] = get_k3D(nx,ny,nz,dx,dy,dz)
Lx = nx*dx; 
Ly = ny*dy;
Lz = nz*dz;
kx = fftshift((-nx/2:nx/2-1)*(2*pi/Lx));
ky = fftshift((-ny/2:ny/2-1)*(2*pi/Ly));
kz = fftshift((-nz/2:nz/2-1)*(2*pi/Lz));
[k1,k2,k3] = ndgrid(kx,ky,kz);
kxdeal = 1/3*nx*(2*pi/Lx);
kydeal = 1/3*ny*(2*pi/Ly);
kzdeal = 1/3*nz*(2*pi/Lz);
