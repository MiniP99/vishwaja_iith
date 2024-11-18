%% Write  self-similar data for Mc-0.2--------
clear,clf
inputdir='/home/vsj/Codes/runs/run02/base/';
outputdir='/home/vsj/Codes/perturb_code/initial02/';
mkdir(outputdir);
nx = 128; ny = 180; nz = 64 ;


time=425;
counter=time/25;
file=strcat(inputdir,'shearlayer1mat_',sprintf('%04d',counter),'.h5');
u=h5read(file,'//u');v=h5read(file,'//v');w=h5read(file,'//w');
rho=h5read(file,'//rho');p=h5read(file,'//p');


Lx = 150; Ly = 200; Lz = 75;
nx = 128; ny = 180; nz = 64 ;
dx = Lx/nx; x = linspace(0,Lx-dx,nx); 
dy = Ly/ny; y = linspace(-Ly/2,Ly/2,ny); 
dz = Lz/nz; z = linspace(0,Lz-dz,nz); 
c=sqrt(1.4);Mc=0.2;du=Mc*2*c;rho_0=1;Re_theta=1000;
p_ref=1;rho_ref=1;
u    = u - ((du/2)*tanh(y/2));
p    = p - p_ref;
rho  = rho - rho_ref;

figure(1), clf
contourf(x, y, squeeze(v(:,:,nz/2))', 'LineStyle', 'none'), pbaspect([Lx Ly 1])
title('u'), xlabel('x'), ylabel('y'), colorbar
screen2jpeg(strcat('vxy_pert.png'))


fname = strcat(outputdir,'perturb_u_',sprintf('%04d',nx),'_',sprintf('%04d',ny),'_',sprintf('%04d',nz),'.dat');
done = write_fortran_box(fname, u, 'double');

fname = strcat(outputdir,'perturb_v_',sprintf('%04d',nx),'_',sprintf('%04d',ny),'_',sprintf('%04d',nz),'.dat');
done = write_fortran_box(fname, v, 'double');

fname = strcat(outputdir,'perturb_w_',sprintf('%04d',nx),'_',sprintf('%04d',ny),'_',sprintf('%04d',nz),'.dat');
done = write_fortran_box(fname, w, 'double');


fname = strcat(outputdir,'perturb_p_',sprintf('%04d',nx),'_',sprintf('%04d',ny),'_',sprintf('%04d',nz),'.dat');
done = write_fortran_box(fname, p, 'double');

fname = strcat(outputdir,'perturb_r_',sprintf('%04d',nx),'_',sprintf('%04d',ny),'_',sprintf('%04d',nz),'.dat');
done = write_fortran_box(fname, rho, 'double');


Lx = 80; Ly = 80; Lz = 40;
nx = 128; ny = 180; nz = 64 ;
dx = Lx/nx; x = linspace(0,Lx-dx,nx); 
dy = Ly/ny; y = linspace(-Ly/2,Ly/2,ny); 
dz = Lz/nz; z = linspace(0,Lz-dz,nz); 
c=sqrt(1.4);Mc=2.0;du=Mc*2*c;rho_0=1;Re_theta=1000;
p_ref=1;rho_ref=1;
u    = u + ((du/2)*tanh(y/2));
p    = p + p_ref;
rho  = rho + rho_ref;
maxu= max(u,[],'all')
figure(2), clf
contourf(x, y, squeeze(v(:,:,nz/2))', 'LineStyle', 'none'), pbaspect([Lx Ly 1])
title('u'), xlabel('x'), ylabel('y'), colorbar
screen2jpeg(strcat('vxy.png'))

