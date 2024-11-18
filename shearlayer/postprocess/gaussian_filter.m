function [f_filz] = gaussian_filter(f,nx,ny,nz)
f_filx= gauss_x(f,nx,ny,nz);
f_fily= gauss_y(f_filx,nx,ny,nz);
f_filz= gauss_z(f_fily,nx,ny,nz);
end
