function done = plot_contours_uvw(x, y, z, nx, ny, nz, Lx, Ly, Lz, u, u_sc, fieldlabel, counter)
 
  done = 0;

  figure(1), clf
  contourf(x, y, squeeze(u(:,:,nz/2))', 'LineStyle', 'none'), pbaspect([Lx Ly 1])
  title('Original'), xlabel('x'), ylabel('y'), colorbar
  screen2jpeg(strcat(fieldlabel,'contxy_orig_',sprintf('%04d',nx),'_',sprintf('%04d',ny),'_',sprintf('%02d',counter),'.png'))
  
  figure(1), clf
  contourf(x, y, squeeze(u_sc(:,:,nz/2))', 'LineStyle', 'none'), pbaspect([Lx Ly 1])
  title('Scaled'), xlabel('x'), ylabel('y'), colorbar
  screen2jpeg(strcat(fieldlabel,'contxy_scal_',sprintf('%04d',nx),'_',sprintf('%04d',ny),'_',sprintf('%02d',counter),'.png'))
  
  figure(1), clf
  contourf(z, y, squeeze(u(nx/2,:,:)), 'LineStyle', 'none'), pbaspect([Lz Ly 1])
  title('Original'), xlabel('z'), ylabel('y'), colorbar
  screen2jpeg(strcat(fieldlabel,'contyz_orig_',sprintf('%04d',ny),'_',sprintf('%04d',nz),'_',sprintf('%02d',counter),'.png'))
  
  figure(1), clf
  contourf(z, y, squeeze(u_sc(nx/2,:,:)), 'LineStyle', 'none'), pbaspect([Lz Ly 1])
  title('Scaled'), xlabel('z'), ylabel('y'), colorbar
  screen2jpeg(strcat(fieldlabel,'contyz_scal_',sprintf('%04d',ny),'_',sprintf('%04d',nz),'_',sprintf('%02d',counter),'.png'))

  figure(1), clf
  contourf(x, z, squeeze(u(:,ny/2,:))', 'LineStyle', 'none'), pbaspect([Lx Lz 1])
  title('Original'), xlabel('x'), ylabel('z'), colorbar
  screen2jpeg(strcat(fieldlabel,'contxz_orig_',sprintf('%04d',nx),'_',sprintf('%04d',nz),'_',sprintf('%02d',counter),'.png'))
  
  figure(1), clf
  contourf(x, z, squeeze(u_sc(:,ny/2,:))', 'LineStyle', 'none'), pbaspect([Lx Lz 1])
  title('Scaled'), xlabel('x'), ylabel('z'), colorbar
  screen2jpeg(strcat(fieldlabel,'contxz_scal_',sprintf('%04d',nx),'_',sprintf('%04d',nz),'_',sprintf('%02d',counter),'.png'))

   done = 1;

end
