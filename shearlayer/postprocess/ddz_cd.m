function [f_p]=ddz_cd(f,dz,nx,ny,nz)

f_p = zeros(nx,ny,nz);
f_p(:,:,1)     = (f(:,:,2) - f(:,:,nz))/(2*dz);
for i=2:nz-1
    f_p(:,:,i) = (f(:,:,i+1) - f(:,:,i-1))/(2*dz);
end
f_p(:,:,nz)     = (f(:,:,1) - f(:,:,nz-1))/(2*dz);
