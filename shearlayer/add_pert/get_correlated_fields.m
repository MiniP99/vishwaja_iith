function u_sc = get_correlated_fields(u, uh, x, y, z, kx, ky, kz, nx, ny, nz, Lx, Ly, Lz, fieldlabel, counter)

  [k3D, Ek3D, ksph, Ek, iksphA, iksphC] = get_1D_spectrum(kx, ky, kz, nx, ny, nz, uh);
  
  figure(1)
  subplot(2,1,1); plot(x, u(:,ny/2,1), '-')
  subplot(2,1,2); loglog(ksph, Ek, '-o')
  screen2jpeg(strcat('sigspec_orig_',sprintf('%04d',nx),'_',sprintf('%02d',counter),'.png'))
  
  % model spectrum
  Ekmod = zeros(size(Ek));
  kPeak =48*2*pi/Lx;
  Ekmod = (ksph/kPeak).^4 .* exp(-2*(ksph/kPeak).^2);
  Ekmod3D = (k3D/kPeak).^4 .* exp(-2*(k3D/kPeak).^2);
  figure(1), clf
  loglog(ksph, Ekmod, '-o')
  screen2jpeg('model_spectrum.png')
  
  % scale the spectrum and find physical field
  uh_scaled = uh.*sqrt(Ekmod3D);
  [k3D_dum, Ek3D_dum, ksph_sc, Ek_sc, iksphA_sc, iksphC_sc] = get_1D_spectrum(kx, ky, kz, nx, ny, nz, uh_scaled);
  u_sc = ifftn(uh_scaled);
  u_sc_check = max(abs(imag(u_sc)./abs(u_sc)),[],'all')
  
  figure(1), clf
  subplot(2,1,1); plot(x, u_sc(:,nz/2,1), '-')
  subplot(2,1,2); loglog(ksph_sc, Ek_sc, '-o')
  screen2jpeg(strcat('sigspec_scal_',sprintf('%04d',nx),'_',sprintf('%02d',counter),'.png'))
 
  done = plot_contours_uvw(x, y, z, nx, ny, nz, Lx, Ly, Lz, u, u_sc, fieldlabel, counter);

end
