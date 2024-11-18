function [f_p]=ddy_cd(f,dy,nx,ny,nz)
%symmetric boundary condition
f_p = zeros(nx,ny,nz);
f_p(:,1,:)     =  0;
for i=2:ny-1
    f_p(:,i,:) = (f(:,i+1,:) - f(:,i-1,:))/(2*dy);
end
f_p(:,ny,:)     = 0;
