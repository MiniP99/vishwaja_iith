function [f_p]=ddx_cd(f,dx,nx,ny,nz)

f_p = zeros(nx,ny,nz);
f_p(1,:,:)     = (f(2,:,:) - f(nx,:,:))/(2*dx);
for i=2:nx-1
    f_p(i,:,:) = (f(i+1,:,:) - f(i-1,:,:))/(2*dx);
end
f_p(nx,:,:)     = (f(1,:,:) - f(nx-1,:,:))/(2*dx);
