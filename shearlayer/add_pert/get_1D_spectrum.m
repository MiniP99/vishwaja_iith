function [k3D, Ek3D, k3Dunique, Ek, ia, ic] = get_1D_spectrum(kx, ky, kz, nx, ny, nz, uh)

  Ek3D = uh.*conj(uh);
  Ek3D = (0.5/(nx*ny*nz))*Ek3D;
  
  
  k3D = zeros(nx, ny, nz);
  for i=1:nx
   for j=1:ny
    for k=1:nz
      k3D(i,j,k) = sqrt(kx(i)^2 + ky(j)^2 + kz(k)^2);
    end
   end
  end
  
  k3D_in_1D = k3D(:);
  Ek3D_in_1D = Ek3D(:);
  
  [k3Dunique, ia, ic] = uniquetol(k3D_in_1D,1e-6);
  nk = length(k3Dunique);
  Ek = zeros(nk,1);
  
  for i=1:length(Ek3D_in_1D)
    ind_to_add = ic(i);
    Ek(ind_to_add) = Ek(ind_to_add) + Ek3D_in_1D(i);
  end

end
