function [u_tilde,u_pp]=favre_avg_fluct(rho,u,nx,ny,nz)
u_tilde=favre(rho,u,nx,nz);
u_pp=fluct(u,u_tilde,nx,ny,nz);
end
