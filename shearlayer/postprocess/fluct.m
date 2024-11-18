function [u_pp]=fluct(u,u_tilde,nx,ny,nz)
u_pp=zeros(nx,ny,nz);
for i=1:nx
    for j=1:ny
        for k=1:nz
            u_pp(i,j,k)=u(i,j,k)-u_tilde(1,j);
        end
    end
end

end
