function [f_filz] = box_filter(f,Nx,Ny,Nz)
f_filx= box_x(f,Nx);
f_fily= box_y(f_filx,Ny);
f_filz= box_z(f_fily,Nz);
end
