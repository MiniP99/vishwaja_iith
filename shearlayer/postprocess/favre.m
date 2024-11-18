function [u_tilde]=favre(rho,u,nx,nz)
u_tilde=avgxz(rho.*u,nx,nz)./avgxz(rho,nx,nz);
end
