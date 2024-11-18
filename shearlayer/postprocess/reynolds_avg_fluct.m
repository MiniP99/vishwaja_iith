function [f_bar,f_reynold_fluct]=reynolds_avg_fluct(f,nx,ny,nz)
f_bar=avgxz(f,nx,nz);
f_reynold_fluct=zeros(nx,ny,nz);
for i=1:nx
    for j=1:ny
        for k=1:nz
            f_reynold_fluct(i,j,k)=f(i,j,k)-f_bar(1,j);
        end
    end
end

end
