clear
outputdir='/home/vsj/Codes/perturb_code/stretch/mc0204/';
mkdir(outputdir);
%For Mc=0.3,1.1  ks=48 in get_correlated_fields
%Lx = 345; Ly = 172; Lz = 86;
%nx = 512; ny = 256; nz = 128;
%For Mc=0.2
Lx = 150;  Ly = 200; Lz = 75;
nx = 180; ny = 240; nz = 84;

dx = Lx/nx; x = linspace(0,Lx-dx,nx); 
dy = Ly/ny; y = linspace(-Ly/2,Ly/2,ny); 
dz = Lz/nz; z = linspace(0,Lz-dz,nz); 

kx = fftshift((-nx/2:nx/2-1)*(2*pi/Lx));
ky = fftshift((-ny/2:ny/2-1)*(2*pi/Ly));
kz = fftshift((-nz/2:nz/2-1)*(2*pi/Lz));


%%%% Step 0 :: Stretched y for now
dx = Lx/nx; xu = linspace(0,Lx-dx,nx); 
dy = Ly/ny; yu = linspace(-Ly/2,Ly/2,ny); 
[X,Y] =meshgrid(yu,xu);
% Streched y
tau = 1;
h = Ly/2;
ystart = -Ly/2;
yfocus = 0;

y_c = yfocus - ystart;
num = 1 + ( (exp(tau) -1) * (y_c/h) );
den = 1 + ( (exp(-tau) -1) * (y_c/h) );
B = (1/(2*tau)) * log(num/den);

dy_str = zeros(1,ny-1);

ystre = zeros(1,ny);
for j=1:ny
y_bar = yu(1,j) - ystart ;
ystre(1,j) = y_c * ( 1 + (sinh(tau * ((y_bar/h)-B)) /sinh(tau*B)) )  + ystart;
end


%%%% Step 1 :: Generate (uncorrelated) random velocity fields
u = randn(nx,ny,nz);   v = randn(nx,ny,nz);    w = randn(nx,ny,nz);


%%%% Step 2 :: Make fields correlated with required spectrum
uh = fftn(u); uc = get_correlated_fields(u, uh, x, ystre, z, kx, ky, kz, nx, ny, nz, Lx, Ly, Lz, 'u', 0); 'Generated correlated u:'
vh = fftn(v); vc = get_correlated_fields(v, vh, x, ystre, z, kx, ky, kz, nx, ny, nz, Lx, Ly, Lz, 'v', 0); 'Generated correlated v:'
wh = fftn(w); wc = get_correlated_fields(w, wh, x, ystre, z, kx, ky, kz, nx, ny, nz, Lx, Ly, Lz, 'w', 0); 'Generated correlated w:'

u = uc; v = vc; w = wc; uc = []; vc = []; wc = [];

% scale [u,v,w] to get the correct TI
TI_reqd = 0.1; U_scale = 0.707;
[u2, v2, w2] = scale_uvw_TI(u, v, w, TI_reqd, U_scale);
u = u2; v = v2; w = w2; u2 = []; v2 = []; w2 = [];

%%%% Step 3 :: Generate pressure field from Poisson equation %%%
[kx3D,ky3D,kz3D,kxdeal,kydeal,kzdeal] = get_k3D(nx,ny,nz,dx,dy,dz);
ksq = (kx3D.^2 + ky3D.^2 + kz3D.^2);

dudx = ddx_compact(u, Lx, nx, ny, nz); dvdy = ddy_compact(v, Ly, nx, ny, nz); dwdz = ddz_compact(w, Lz, nx, ny, nz);
dilat = dudx + dvdy + dwdz;
dilath = fftn(dilat); phat = -dilath./ksq; phat(1,1,1) = 0; p = ifftn(phat);
p_check = max(abs(imag(p)./abs(p)),[],'all')

%%%% Step 4 :: Make velocity field solenoidal %%%
dpdx = ddx(p, kx3D);   dpdy = ddy(p, ky3D);   dpdz = ddz(p, kz3D);
usol = u - dpdx; vsol = v - dpdy; wsol = w - dpdz;

maxdiv = get_div(u, v, w, kx3D, ky3D, kz3D);  % check divergence of [u, v, w]
maxdivsol = get_div(usol, vsol, wsol, kx3D, ky3D, kz3D);  % check divergence of [usol, vsol, wsol]
[maxdiv maxdivsol]

% fix TI level
TIsol_current = sqrt(mean(usol.^2 + vsol.^2 + wsol.^2,'all')/3)/U_scale; TIsol_current
[u2, v2, w2] = scale_uvw_TI(usol, vsol, wsol, TI_reqd, U_scale);
usol = u2; vsol = v2; wsol = w2; u2 = []; v2 = []; w2 = [];

% check solenoidality again
maxdivsol = get_div(usol, vsol, wsol, kx3D, ky3D, kz3D);  % check divergence of [usol, vsol, wsol]
[maxdiv maxdivsol]
TIsol_current = sqrt(mean(usol.^2 + vsol.^2 + wsol.^2,'all')/3)/U_scale; TIsol_current

%%%% Step 5 :: Generate density field from pressure field assuming isentropic fluctuations %%%
p_base = 1; rho_base = 1; gam = 1.4;
const = p_base/rho_base^gam;
p = p + p_base;
psol=p-p_base;
rho = (p/const).^(1/gam); %given
rhosol = rho-rho_base;


%%%% Step 6 :: Write all fields to file %%%
%done = plot_contours_uvw(x, y, z, nx, ny, nz, Lx, Ly, Lz, u, usol, 'usol', 0);
%done = plot_contours_uvw(x, y, z, nx, ny, nz, Lx, Ly, Lz, v, vsol, 'vsol', 0);
%done = plot_contours_uvw(x, y, z, nx, ny, nz, Lx, Ly, Lz, w, wsol, 'wsol', 0);
%done = plot_contours_uvw(x, y, z, nx, ny, nz, Lx, Ly, Lz, p, psol,    'psol', 0);
%done = plot_contours_uvw(x, y, z, nx, ny, nz, Lx, Ly, Lz, rho, rhosol,'rhosol', 0);

fname = strcat(outputdir,'perturb_u_',sprintf('%04d',nx),'_',sprintf('%04d',ny),'_',sprintf('%04d',nz),'.dat');
done = write_fortran_box(fname, usol, 'double');

fname = strcat(outputdir,'perturb_v_',sprintf('%04d',nx),'_',sprintf('%04d',ny),'_',sprintf('%04d',nz),'.dat');
done = write_fortran_box(fname, vsol, 'double');

fname = strcat(outputdir,'perturb_w_',sprintf('%04d',nx),'_',sprintf('%04d',ny),'_',sprintf('%04d',nz),'.dat');
done = write_fortran_box(fname, wsol, 'double');


fname = strcat(outputdir,'perturb_p_',sprintf('%04d',nx),'_',sprintf('%04d',ny),'_',sprintf('%04d',nz),'.dat');
done = write_fortran_box(fname, psol, 'double');

fname = strcat(outputdir,'perturb_r_',sprintf('%04d',nx),'_',sprintf('%04d',ny),'_',sprintf('%04d',nz),'.dat');
done = write_fortran_box(fname, rhosol, 'double');
