function [f_bar]=avgxz(f,nx,nz)
avgx=sum(f,1)/nx;
f_bar=sum(avgx,3)/nz;
end
