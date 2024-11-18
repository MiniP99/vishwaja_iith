function [f_bar]=avgyz(f,ny,nz)
avgx=sum(f,2)/ny;
f_bar=sum(avgx,3)/nz;
end
