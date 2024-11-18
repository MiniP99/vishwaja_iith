function [R11,R12,R13,R22,R23,R33]=reynold_stress(rho,u_pp,v_pp,w_pp,nx,nz)
R11=favre(rho,u_pp.*u_pp,nx,nz);
R12=favre(rho,u_pp.*v_pp,nx,nz);
R13=favre(rho,u_pp.*w_pp,nx,nz);
R22=favre(rho,v_pp.*v_pp,nx,nz);
R23=favre(rho,v_pp.*w_pp,nx,nz);
R33=favre(rho,w_pp.*w_pp,nx,nz);
end

